import json
import sys
from collections import defaultdict
import pandas as pd
from Bio import SeqIO
from functional import pseq
from src.utils import *


def main(dest_dir, threads=4, ident_min=95):
    print("## 01_reformat_ffn ##")

    ffn_dir = joinpath(source_dir, "FFN")
    ffn2_dir = joinpath(dest_dir, "FFN2")
    create_if_not_exist(ffn2_dir)

    ffn_files = os.listdir(ffn_dir)
    for filename in ffn_files:
        handle = SeqIO.parse(joinpath(ffn_dir, filename), "fasta")
        records = [">" + record.id + "\n" + str(record.seq) for record in handle]

        with open(joinpath(ffn2_dir, filename), "w") as file:
            file.write("\n".join(records))

    print("--Finished")


    # calculate the pan genome
    print("## 02_runRoary ##")

    cmd1 = "docker run --rm -v {dest_dir}:/data -v {source_dir}/GFF:/GFF a504082002/roary".format(**locals())
    cmd2 = "roary -e --mafft -p {threads} -i {ident_min} -f /data/roary /GFF/*.gff".format(**locals())
    os.system(cmd1 + " " + cmd2)

    # change owner
    shutil.chown(joinpath(dest_dir, "roary"), user=os.getuid())

    print("--Finished")



    print("## 03_parse_csvfile ##")

    roary_dir = joinpath(dest_dir, "roary")
    locus_dir = joinpath(dest_dir, "locusfiles")
    create_if_not_exist(locus_dir)

    roary_result = pd.read_csv(joinpath(roary_dir, "gene_presence_absence.csv"))
    roary_result.index = pd.Index(map(lambda x: "SAL{0:07d}".format(x+1), roary_result.index), name="locus")

    mapping = roary_result[roary_result["No. isolates"] == roary_result["No. sequences"]]
    mapping = mapping[["Gene", "No. isolates", "No. sequences", "Annotation"]]
    mapping.to_csv(joinpath(dest_dir, "locusmapping.txt"), sep="\t")

    paralog = roary_result[roary_result["No. isolates"] != roary_result["No. sequences"]]
    paralog = paralog[["Gene", "No. isolates", "No. sequences", "Annotation"]]
    paralog.to_csv(joinpath(dest_dir, "paralog.txt"), sep="\t", index=False)

    for colname, col in roary_result.iloc[:, 14:].iteritems():
        handle = SeqIO.parse(joinpath(ffn2_dir, colname + ".ffn"), "fasta")
        cds = {record.id: record.seq for record in handle}

        for rowname, loci in col.iteritems():
            if isinstance(loci, str):
                with open(joinpath(locus_dir, rowname + ".tmp"), "a") as file:
                    prokkaID = loci.split("___")[0]
                    file.write(">" + prokkaID + "\n" + str(cds[prokkaID]) + "\n")

    print("--Finished")

    print("## 04_collapse ##")

    locus_files = os.listdir(locus_dir)
    pseq(locus_files).map(run_fastx)

    # change owner
    os.system("chown -R $USER:users " + locus_dir)

    print("--Finished")


    print("## 05_renumberAlleles ##")

    for locus_file in locus_files:
        filename = locus_file.split(".")[0]
        file_path = joinpath(dest_dir, "locusfiles", filename)

        records = []
        for record in SeqIO.parse(file_path, "fasta"):
            new_id = record.id.split("-")[0]
            new_seq = str(record.seq)
            records.append(">{new_id}\n{new_seq}".format(**locals()))

        with open(file_path, "w") as file:
            file.write("\n".join(records))

    print("--Finished")


    print("## 06_makePanRef ##")

    # find max length sequence in a locus file as the pan RefSeq
    refseqs = {}
    for locus_file in locus_files:
        filename = locus_file.split(".")[0]
        maxlen_seq = ""
        maxlen = 0
        for record in SeqIO.parse(joinpath(dest_dir, "locusfiles", filename), "fasta"):
            current_seq = str(record.seq)
            current_len = len(current_seq)
            if current_len > maxlen:
                maxlen = current_len
                maxlen_seq = current_seq

        refseqs[filename] = maxlen_seq

    with open(joinpath(dest_dir, "panRefSeq.fa"), "w") as file:
        seqs = [">{key}\n{value}".format(**locals()) for key, value in refseqs.items()]
        file.write("\n".join(seqs))

    with open(joinpath(dest_dir, "panRefSeq.json"), "w") as file:
        file.write(json.dumps(refseqs))

    print("--Finished")



    # 07_makeSchems

    print("## 07_makeSchems ##")

    scheme_dir = dest_dir + "/scheme"
    if not os.path.exists(scheme_dir):
        os.mkdir(scheme_dir)

    totisolates = len(os.listdir(ffn2_dir))

    dispatcher = defaultdict(list)
    for index, row in pd.read_csv(dest_dir + "/locusmapping.txt", sep="\t").iterrows():
        locus = row["locus"]
        occ = row["No. isolates"]
        occrate = occ / totisolates

        dispatcher["pan"].append(locus)
        if 0.95 <= occrate <= 1:
            dispatcher["core"].append(locus)
        if 0.85 <= occrate <= 1:
            dispatcher["core_10p_dispensable"].append(locus)
        if 0.75 <= occrate <= 1:
            dispatcher["core_20p_dispensable"].append(locus)
        if 0.65 <= occrate <= 1:
            dispatcher["core_30p_dispensable"].append(locus)
        if 0.55 <= occrate <= 1:
            dispatcher["core_40p_dispensable"].append(locus)
        if 0.45 <= occrate <= 1:
            dispatcher["core_50p_dispensable"].append(locus)
        if 0.35 <= occrate <= 1:
            dispatcher["core_60p_dispensable"].append(locus)
        if 0.25 <= occrate <= 1:
            dispatcher["core_70p_dispensable"].append(locus)
        if 0.15 <= occrate <= 1:
            dispatcher["core_80p_dispensable"].append(locus)
        if 0.05 <= occrate <= 1:
            dispatcher["core_90p_dispensable"].append(locus)
        if occ > 1 and occrate <= 1:
            dispatcher["core_100p_dispensable"].append(locus)
        if occ > 1 and occrate <= 0.95:
            dispatcher["dispensable"].append(locus)
        if occ == 1:
            dispatcher["unique"].append(locus)

    with open(scheme_dir + "/scheme.txt", "w") as file2:
        for key, loci in dispatcher.items():
            if len(loci) != 0:
                with open(scheme_dir + "/" + key + ".txt", "w") as file:
                    content = [">" + l + "\n" + refseqs[l] for l in loci]
                    file.write("\n".join(content))
                content2 = key + "\t" + "\t".join(loci) + "\n"
                file2.write(content2)

    print("--Finished")



    # 08_output

    print("## 08_output ##")

    output_dir = dest_dir + "/OUTPUT"
    os.mkdir(output_dir)
    shutil.copy(dest_dir + "/roary/summary_statistics.txt", output_dir)
    shutil.copy(dest_dir + "/locusmapping.txt", output_dir)
    shutil.move(scheme_dir, output_dir)
    shutil.copytree(locus_dir, dest_dir + "/locusfiles.bak")
    for file in locus_files:
        if file.endswith(".tmp"):
            os.remove(os.path.join(locus_dir, file))
    shutil.move(locus_dir, output_dir)

    print("--Finished")


    for filename in os.listdir(output_dir + "/locusfiles"):
        handle = SeqIO.parse(output_dir + "/locusfiles/" + filename, "fasta")
        with open(output_dir + "/locusfiles/" + filename + ".new", "w") as file:
            records = [">" + filename + "::" + record.id + "\n" + str(record.seq) + "\n" for record in handle]
            file.write("".join(records))


    print("Congradulation!! All jobs are finished!!")


if __name__ == "__main__":
    root_dir = "/media/pika/Workbench/workspace/forYH_MSGAdb-builder"
    dest_dir = root_dir + "/02--BuidProfileDB/test"
    source_dir = root_dir + "/01--ContigAnnotation/testInputFiles"
    dest_dir = sys.argv[1]
    main(dest_dir)
