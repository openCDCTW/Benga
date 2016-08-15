import sys
import json
from collections import defaultdict, OrderedDict
from functional import seq, pseq

import pandas as pd
import matplotlib.pyplot as plt
import fastcluster
from scipy.cluster.hierarchy import dendrogram
from Bio import SeqIO

from utils import *


def annotate_configs(dest_dir, source_dir):
    print("## 02_reformatContigFile ##")

    filenames = parse_filenames(source_dir)
    working_dir = joinpath(dest_dir, "Assembly")
    create_if_not_exist(working_dir)

    namemap = {}
    for i, oldname in enumerate(filenames, 1):
        newname = "Assembly_" + str(i) + ".fa"
        namemap[oldname] = newname
        with open(joinpath(working_dir, newname), "w") as file:
            for j, contig in enumerate(SeqIO.parse(joinpath(source_dir, oldname), "fasta"), 1):
                seqid = "A_{i}::C_{j}".format(**locals())
                SeqIO.write(replace_id(contig, seqid), file, "fasta")

    print("--Finished")


    print("## 03_contigAnn ##")

    annotate_dir = joinpath(dest_dir, "Assembly_ann")
    create_if_not_exist(annotate_dir)
    pseq(namemap.values()).map(lambda x: annotate(x, dest_dir)).to_list()

    # change owner
    os.chown(annotate_dir, os.getuid(), os.getgid())

    print("--Finished")


    print("## 04_make_ffn-gff ##")

    gff_dir = joinpath(dest_dir, "GFF")
    ffn_dir = joinpath(dest_dir, "FFN")
    create_if_not_exist(gff_dir)
    create_if_not_exist(ffn_dir)

    for folder in os.listdir(annotate_dir):
        path = joinpath(annotate_dir, folder)
        for file in os.listdir(path):
            current_file = joinpath(path, file)
            if file.endswith(".ffn"):
                os.rename(current_file, joinpath(ffn_dir, file))
            elif file.endswith(".gff"):
                os.rename(current_file, joinpath(gff_dir, file))

    print("--Finished")


    print("## 05_storeNonCDS ##")

    noncds = defaultdict(list)
    for file in os.listdir(gff_dir):
        name, ext = file.split(".")
        for line in open(joinpath(gff_dir, file), "r").read().splitlines():
            if line.startswith("##sequence") or line.startswith("##gff"):
                continue
            if line.startswith("##FASTA"):
                break

            token = line.split("\t")
            seq_type, annotation = token[2], token[8]

            if annotation.startswith("ID="):
                prokkaid = annotation.split(";")[0][3:]
                if seq_type != "CDS":
                    noncds[name].append(prokkaid)

    with open(joinpath(dest_dir, "nonCDS.json"), "w") as file:
        file.write(json.dumps(noncds))

    print("--Finished")


def make_profiles(dest_dir, threads=2, ident_min=95):
    print("## 01_reformat_ffn ##")

    ffn_dir = joinpath(dest_dir, "FFN")
    ffn2_dir = joinpath(dest_dir, "FFN2")
    create_if_not_exist(ffn2_dir)

    for filename in os.listdir(ffn_dir):
        handle = SeqIO.parse(joinpath(ffn_dir, filename), "fasta")
        records = [new_record(record.id, str(record.seq)) for record in handle]
        SeqIO.write(records, joinpath(ffn2_dir, filename), "fasta")

    print("--Finished")


    # calculate the pan genome
    print("## 02_runRoary ##")

    cmd1 = "docker run --rm -v {dest_dir}:/data a504082002/roary".format(**locals())
    cmd2 = "python /program/cmd.py roary -e --mafft -p {threads} -i {ident_min} -f /data/roary /data/GFF/*.gff".format(**locals())
    os.system(cmd1 + " " + cmd2)

    # change owner
    os.chown(joinpath(dest_dir, "roary"), os.getuid(), os.getgid())

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
    files_per_docker = 500
    pseq(
        partition(locus_files, files_per_docker)
    ).map(
        lambda x: run_fastx(x, locus_dir)
    ).to_list()

    # change owner
    os.chown(locus_dir, os.getuid(), os.getgid())

    print("--Finished")


    print("## 05_renumberAlleles ##")

    for locus_file in locus_files:
        file_path = joinpath(dest_dir, "locusfiles", locus_file.split(".")[0])

        records = seq(
            SeqIO.parse(file_path, "fasta")
        ).map(
            lambda x: new_record(x.id.split("-")[0], str(x.seq))
        ).to_list()

        SeqIO.write(records, file_path, "fasta")

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

    records = [new_record(key, value) for key, value in refseqs.items()]
    SeqIO.write(records, joinpath(dest_dir, "panRefSeq.fa"), "fasta")

    with open(joinpath(dest_dir, "panRefSeq.json"), "w") as file:
        file.write(json.dumps(refseqs))

    print("--Finished")



    # 07_makeSchems

    print("## 07_makeSchems ##")

    scheme_dir = joinpath(dest_dir, "scheme")
    create_if_not_exist(scheme_dir)

    totisolates = len(os.listdir(ffn2_dir))
    dispatcher = defaultdict(list)
    for index, row in pd.read_csv(joinpath(dest_dir, "locusmapping.txt"), sep="\t").iterrows():
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

    schemes = []
    for key, loci in dispatcher.items():
        if len(loci) != 0:
            records = [new_record(l, refseqs[l]) for l in loci]
            SeqIO.write(records, joinpath(scheme_dir, key + ".txt"), "fasta")
            schemes.append(key + "\t" + "\t".join(loci))

    with open(joinpath(scheme_dir, "scheme.txt"), "w") as file2:
        file2.write("\n".join(schemes))

    print("--Finished")



    # 08_output

    print("## 08_output ##")

    output_dir = joinpath(dest_dir, "OUTPUT")
    os.mkdir(output_dir)
    shutil.copy(joinpath(dest_dir, "roary", "summary_statistics.txt"), output_dir)
    shutil.copy(joinpath(dest_dir, "locusmapping.txt"), output_dir)
    shutil.move(scheme_dir, output_dir)
    shutil.copytree(locus_dir, joinpath(dest_dir, "locusfiles.bak"))
    for file in locus_files:
        if file.endswith(".tmp"):
            os.remove(joinpath(locus_dir, file))
    shutil.move(locus_dir, output_dir)

    print("--Finished")

    locusfiles_dir = joinpath(output_dir, "locusfiles")
    for filename in os.listdir(locusfiles_dir):
        handle = SeqIO.parse(joinpath(locusfiles_dir, filename), "fasta")
        records = [new_record(filename + "::" + record.id, str(record.seq)) for record in handle]
        SeqIO.write(records, joinpath(locusfiles_dir, filename + ".new"), "fasta")

    print("Congratulation!! All jobs are finished!!")


def wg_profiling(dest_dir, blastn_threads=2, aligcov_cut=0.5, pident_cut=90):
    # 01_reformatContigFile
    print("## 01_reformatContigFile ##")

    assemble_dir = joinpath(dest_dir, "Assembly")
    input_dir = joinpath(dest_dir, "input")
    create_if_not_exist(assemble_dir)

    namemap = {}
    for i, filename in enumerate(sorted(os.listdir(input_dir)), 1):
        records = []
        file = SeqIO.parse(joinpath(input_dir, filename), "fasta")
        for j, record in enumerate(file, 1):
            newid = "Assembly_{i}::Contig_{j}".format(**locals())
            new_record = SeqRecord(Seq(str(record.seq), generic_dna), id=newid, description="")
            records.append(new_record)

        newname = "Assembly_{i}.fa".format(**locals())
        SeqIO.write(records, joinpath(assemble_dir, newname), "fasta")
        namemap[newname] = filename

    print("--Finished")

    # 02_locusProfiling
    print("## 02_locusProfiling ##")

    scheme_dir = joinpath(dest_dir, "scheme")
    blastquery = joinpath(scheme_dir, "pan.txt")
    shutil.copytree(joinpath(source_dir, "scheme"), scheme_dir)

    seqlen = {record.id: len(record.seq) for record in SeqIO.parse(blastquery, "fasta")}
    m = [(recordid, "") for recordid in seqlen.keys()]
    m = OrderedDict(sorted(m, key=lambda x: x[0]))

    blast_names = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend",
                   "sstart", "send", "evalue", "bitscore"]
    blast_items = " ".join(blast_names)
    for filename in os.listdir(assemble_dir):
        os.system("makeblastdb -in {assemble_dir}/{filename} -dbtype nucl -out {dest_dir}/lociDB".format(**locals()))
        os.system(
            "blastn -num_threads {blastn_threads} -query {blastquery} -db {dest_dir}/lociDB -out {dest_dir}/blast.out -outfmt '6 {blast_items}'".format(
                **locals()))

        basename = filename.split(".")[0]
        c = {record.id: str(record.seq) for record in SeqIO.parse(joinpath(assemble_dir, filename), "fasta")}
        clone_m = m.copy()

        blast = pd.read_csv(joinpath(dest_dir, "blast.out"), sep="\t", header=None, names=blast_names)
        for index, row in blast.iterrows():
            qid = row["qseqid"]
            qlen = seqlen[qid]
            aligcov = (row["length"] - row["gapopen"]) / qlen
            if aligcov >= aligcov_cut and row["pident"] >= pident_cut:
                row_elements = [aligcov, row["length"], qlen, row["pident"], row["sseqid"], row["sstart"], row["send"]]
                clone_m[qid] = "\t".join(map(str, row_elements))

        newname = "locusAP.{basename}.list".format(**locals())
        with open(joinpath(dest_dir, newname), "w") as file:
            for key, value in clone_m.items():
                hit = 1 if len(value) > 1 else 0
                file.write("{key}\t{hit}\n".format(**locals()))

        os.remove(joinpath(dest_dir, "blast.out"))

    print("--Finished")

    print("## 03_alleleAllocate ##")

    for filename in os.listdir(assemble_dir):
        qfn2 = filename.split(".")[0]
        print(qfn2)
        os.mkdir(joinpath(dest_dir, qfn2))

        blastinfile = joinpath(assemble_dir, filename)
        os.system("makeblastdb -in {blastinfile} -dbtype nucl -out {dest_dir}/AssemblyDB_{qfn2}".format(**locals()))
        assdic = {record.id: str(record.seq) for record in SeqIO.parse(blastinfile, "fasta")}

        allprofile = defaultdict(list)
        newloc = []
        aldic = {}
        for line in open("{dest_dir}/locusAP.{qfn2}.list".format(**locals()), "r").read().splitlines():
            locus, hits = line.split("\t")
            if hits == "0":
                allprofile[locus] = [0]
                newloc.append(locus)
            else:
                qfile = "{source_dir}/locusfiles/{locus}.new".format(**locals())
                for record in SeqIO.parse(qfile, "fasta"):
                    alno = record.id.split("::")[1]
                    aldic[(locus, alno)] = len(str(record.seq))
                os.system("cat {qfile} >> {dest_dir}/qfile.fa".format(**locals()))

        os.system(
            "blastn -num_threads {blastn_threads} -query {dest_dir}/qfile.fa -db {dest_dir}/AssemblyDB_{qfn2} -out {dest_dir}/blast.{qfn2}.out -outfmt '6 {blast_items}'".format(
                **locals()))
        blast = pd.read_csv("{dest_dir}/blast.{qfn2}.out".format(**locals()), sep="\t", header=None, names=blast_names)
        blast["loc"] = blast["qseqid"].str.split("::").str[0]
        blast["allno"] = blast["qseqid"].str.split("::").str[1]

        for _, row in blast[(blast["pident"] == 100) & (blast["mismatch"] == 0) & (blast["gapopen"] == 0)].iterrows():
            if row["length"] == aldic[(row["loc"], row["allno"])]:
                allprofile[row["loc"]].append((row["allno"], row["length"]))

        # New alleles
        new_alleles = []
        records = []
        for _, row in blast[~blast["loc"].isin(allprofile.keys())].iterrows():
            allprofile[row["loc"]] = [0]
            new_alleles.append(row["loc"])
            records.append(new_record(row["loc"], assdic[row["sseqid"]]))

        SeqIO.write(drop_duplicate(records, lambda x: x.id), joinpath(dest_dir, qfn2, "tot.1.new"), "fasta")
        with open(joinpath(dest_dir, qfn2, "tot.new"), "w") as file:
            file.write("\n".join(drop_duplicate(new_alleles)))

        output = []
        for locus, value in allprofile.items():
            if all(map(lambda x: x == 0, value)):
                a = (locus, 0)
            else:
                # only access the last of max length allele number in values
                a = (locus, max(reversed(value), key=lambda x: x[1])[0])
            output.append(a)

        with open(joinpath(dest_dir, "allele." + qfn2 + ".profile"), "w") as file:
            lines = list(map(lambda x: str(x[0]) + "\t" + str(x[1]), sorted(output, key=lambda x: x[0])))
            file.write("\n".join(lines))
        os.remove(joinpath(dest_dir, "qfile.fa"))

    print("--Finished")

    print("## 04_catchSchemeProfile ##")

    scheme = {}
    for line in open(joinpath(scheme_dir, "scheme.txt"), "r").read().splitlines():
        token = line.split("\t")
        scheme[token[0]] = token[1:]
        os.mkdir(joinpath(dest_dir, "locusAP_" + token[0]))
        os.mkdir(joinpath(dest_dir, "wgMLST_" + token[0]))

    for filename in os.listdir(assemble_dir):
        ass = filename.split(".")[0]
        for key, value in scheme.items():
            locusAP_output = seq(
                open(joinpath(dest_dir, "locusAP." + ass + ".list"), "r").read().splitlines()
            ).map(
                lambda x: x.split("\t")
            ).filter(
                lambda x: x[0] in value
            ).map(
                lambda x: x[0] + "\t" + x[1]
            ).to_list()

            wgMLST_output = seq(
                open(joinpath(dest_dir, "allele." + ass + ".profile"), "r").read().splitlines()
            ).map(
                lambda x: x.split("\t")
            ).filter(
                lambda x: x[0] in value
            ).map(
                lambda x: x[0] + "\t" + x[1]
            ).to_list()

            with open("{dest_dir}/locusAP_{key}/locusAP.{ass}.{key}.list".format(**locals()), "w") as file:
                file.write("\n".join(locusAP_output))
            with open("{dest_dir}/wgMLST_{key}/allele.{ass}.{key}.list".format(**locals()), "w") as file:
                file.write("\n".join(wgMLST_output))

    print("--Finished")

    # 05_outformat
    print("## 05_outformat ##")

    output_dir = joinpath(dest_dir, "OUTPUT")
    clear_folder(output_dir)

    locusap_dir = joinpath(output_dir, "locusAP")
    clear_folder(locusap_dir)

    wgmlst_dir = joinpath(output_dir, "wgMLST")
    clear_folder(wgmlst_dir)

    os.system("mv {dest_dir}/locusAP_* {locusap_dir}".format(**locals()))
    os.system("mv {dest_dir}/wgMLST_* {wgmlst_dir}".format(**locals()))

    dir_names = ["_core", "_pan", "_core_10p_dispensable", "_core_20p_dispensable", "_core_30p_dispensable",
                 "_core_40p_dispensable", "_core_50p_dispensable", "_core_60p_dispensable", "_core_70p_dispensable",
                 "_core_80p_dispensable", "_core_90p_dispensable", "_core_100p_dispensable"]
    for name in dir_names:
        locusap_name = joinpath(locusap_dir, "locusAP" + name)
        os.system("ls {locusap_name} | grep 'locusAP'> {locusap_name}/files.list".format(**locals()))
        wgmlst_name = joinpath(wgmlst_dir, "wgMLST" + name)
        os.system("ls {wgmlst_name} | grep 'allele'> {wgmlst_name}/files.list".format(**locals()))

    os.system("rm -rf Assembly*")

    print("--Finished")


def build_tree(dest_dir, scheme_selection="locusAP_core_100p_dispensable"):
    # cpScheme
    # shutil.copy(joinpath(src_dir, "files.new.list"), ".")
    # shutil.copytree(joinpath(src_dir, "OUTPUT"), "scheme")

    # makeDiffMatrix

    scheme_path = {"locusAP_core": "/scheme/locusAP/locusAP_core",
                   "locusAP_core_10p_dispensable": "/scheme/locusAP/locusAP_core_10p_dispensable",
                   "locusAP_core_20p_dispensable": "/scheme/locusAP/locusAP_core_20p_dispensable",
                   "locusAP_core_30p_dispensable": "/scheme/locusAP/locusAP_core_30p_dispensable",
                   "locusAP_core_40p_dispensable": "/scheme/locusAP/locusAP_core_40p_dispensable",
                   "locusAP_core_50p_dispensable": "/scheme/locusAP/locusAP_core_50p_dispensable",
                   "locusAP_core_60p_dispensable": "/scheme/locusAP/locusAP_core_60p_dispensable",
                   "locusAP_core_70p_dispensable": "/scheme/locusAP/locusAP_core_70p_dispensable",
                   "locusAP_core_80p_dispensable": "/scheme/locusAP/locusAP_core_80p_dispensable",
                   "locusAP_core_90p_dispensable": "/scheme/locusAP/locusAP_core_90p_dispensable",
                   "locusAP_core_100p_dispensable": "/scheme/locusAP/locusAP_core_100p_dispensable",
                   "locusAP_pan": "/scheme/locusAP/locusAP_pan",
                   "locusAP_unique": "/scheme/locusAP/locusAP_unique",
                   "wgMLST_core": "/scheme/wgMLST/wgMLST_core",
                   "wgMLST_core_10p_dispensable": "/scheme/wgMLST/wgMLST_core_10p_dispensable",
                   "wgMLST_core_20p_dispensable": "/scheme/wgMLST/wgMLST_core_20p_dispensable",
                   "wgMLST_core_30p_dispensable": "/scheme/wgMLST/wgMLST_core_30p_dispensable",
                   "wgMLST_core_40p_dispensable": "/scheme/wgMLST/wgMLST_core_40p_dispensable",
                   "wgMLST_core_50p_dispensable": "/scheme/wgMLST/wgMLST_core_50p_dispensable",
                   "wgMLST_core_60p_dispensable": "/scheme/wgMLST/wgMLST_core_60p_dispensable",
                   "wgMLST_core_70p_dispensable": "/scheme/wgMLST/wgMLST_core_70p_dispensable",
                   "wgMLST_core_80p_dispensable": "/scheme/wgMLST/wgMLST_core_80p_dispensable",
                   "wgMLST_core_90p_dispensable": "/scheme/wgMLST/wgMLST_core_90p_dispensable",
                   "wgMLST_core_100p_dispensable": "/scheme/wgMLST/wgMLST_core_100p_dispensable",
                   "wgMLST_pan": "/scheme/wgMLST/wgMLST_pan",
                   "wgMLST_unique": "/scheme/wgMLST/wgMLST_unique"}

    target_dir = dest_dir + scheme_path[scheme_selection]

    namemap = {}
    for file in open("files.new.list", "r").read().splitlines():
        a1, a2 = file.split("\t")[0:2]
        newname = a1.split(".")[0]
        oldname = a2.split(".")[0]
        namemap[newname] = oldname

    # collect labels
    nodeLabel = []
    for filename in open(joinpath(target_dir, "files.list"), "r").read().splitlines():
        fn = filename.split(".")[1]
        nodeLabel.append(fn.replace("Assembly_", "A") + "\t" + namemap[fn])

    with open("nodeLableMatch.txt", "w") as file:
        file.write("\n".join(nodeLabel))

    # collect profiles
    matrix = pd.DataFrame()
    for filename in open(joinpath(target_dir, "files.list"), "r").read().splitlines():
        column = pd.Series()
        for line in open(target_dir + "/" + filename).read().splitlines():
            if line == "":
                continue
            locus, hit = line.split("\t")
            column[locus] = True if hit != 0 else False
        sample = filename.split(".")[1]
        matrix[sample] = column

    matrix.to_csv("profiles.tsv", sep="\t")

    Y = fastcluster.linkage(matrix, method="average", metric="hamming")
    dendrogram(Y)
    plt.show()
    plt.close()


def main(dest_dir, source_dir):
    annotate_configs(dest_dir, source_dir)
    make_profiles(dest_dir)


if __name__ == "__main__":
    # dest_dir = "/media/pika/Workbench/workspace/forYH_MSGAdb-builder/01--ContigAnnotation/testInputFiles"
    # source_dir = joinpath(dest_dir, "upload_Assembly")
    dest_dir = sys.argv[1]
    source_dir = sys.argv[2]
    main(dest_dir, source_dir)

