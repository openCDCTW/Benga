import shutil
from collections import OrderedDict, defaultdict
import pandas as pd
from Bio import SeqIO
from src.utils import *

root_dir = "/media/pika/Workbench/workspace/forYH_MSGAdb-builder"
dest_dir = root_dir + "/03--wgProfileing/test"
source_dir = root_dir + "/02--BuidProfileDB/test/OUTPUT"
blastn_threads = 2
aligcov_cut = 0.5
pident_cut = 90


# 01_reformatContigFile
print("## 01_reformatContigFile ##")

assemble_dir = joinpath(dest_dir, "Assembly")
create_if_not_exist(assemble_dir)

input_dir = joinpath(dest_dir, "input")
for i, filename in enumerate(sorted(os.listdir(input_dir)), 1):
    records = []
    file = SeqIO.parse(joinpath(input_dir, filename), "fasta")
    for j, record in enumerate(file, 1):
        newid = "Assembly_{i}::Contig_{j}".format(**locals())
        new_record = SeqRecord(Seq(str(record.seq), generic_dna), id=newid, description="")
        records.append(new_record)

    newname = "Assembly_{i}.fa".format(**locals())
    SeqIO.write(records, joinpath(assemble_dir, newname), "fasta")

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
    os.system("blastn -num_threads {blastn_threads} -query {blastquery} -db {dest_dir}/lociDB -out {dest_dir}/blast.out -outfmt '6 {blast_items}'".format(**locals()))

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

    os.system("blastn -num_threads {blastn_threads} -query {dest_dir}/qfile.fa -db {dest_dir}/AssemblyDB_{qfn2} -out {dest_dir}/blast.{qfn2}.out -outfmt '6 {blast_items}'".format(**locals()))
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
        record = SeqRecord(Seq(assdic[row["sseqid"]], generic_dna), id=row["loc"], description="")
        records.append(record)

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
        output = []
        for line in open("{dest_dir}/locusAP.{ass}.list".format(**locals()), "r").read().splitlines():
            token = line.split("\t")
            if token[0] in value:
                output.append(token[0] + "\t" + token[1])
        with open("{dest_dir}/locusAP_{key}/locusAP.{ass}.{key}.list".format(**locals()), "w") as file:
            file.write("\n".join(output))
        
        for line in open("{dest_dir}/allele.{ass}.profile".format(**locals()), "r").read().splitlines():
            token = line.split("\t")
            if token[0] in value:
                output.append(token[0] + "\t" + token[1])
        with open("{dest_dir}/wgMLST_{key}/allele.{ass}.{key}.list".format(**locals()), "w") as file:
            file.write("\n".join(output))

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

