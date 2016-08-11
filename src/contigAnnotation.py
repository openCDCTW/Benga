import sys
import json
from collections import defaultdict
from Bio import SeqIO
from functional import seq, pseq
from utils import *


def main(dest_dir, source_dir):
    print("## 01_catchInputFileName ##")

    filenames = parse_filenames(source_dir)

    print("--Finished")


    print("## 02_reformatContigFile ##")

    working_dir = joinpath(dest_dir, "Assembly")
    create_if_not_exist(working_dir)

    namemap = {}
    for i, oldname in enumerate(filenames, 1):
        newname = "Assembly_" + str(i) + ".fa"
        namemap[oldname] = newname
        with open(joinpath(working_dir, newname), "w") as file:
            for j, contig in enumerate(SeqIO.parse(os.path.join(source_dir, oldname), "fasta"), 1):
                seqid = "A_{i}::C_{j}".format(**locals())
                SeqIO.write(replace_id(contig, seqid), file, "fasta")

    print("--Finished")


    print("## 03_contigAnn ##")

    annotate_dir = joinpath(dest_dir, "Assembly_ann")
    create_if_not_exist(annotate_dir)
    pseq(namemap.values()).map(lambda x: annotate(x, dest_dir)).to_list()

    # change owner
    os.chown(annotate_dir, os.getuid(), os.getgid())
    # os.system("chown -R $USER:users " + annotate_dir)

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

if __name__ == "__main__":
    # dest_dir = "/media/pika/Workbench/workspace/forYH_MSGAdb-builder/01--ContigAnnotation/testInputFiles"
    # source_dir = joinpath(dest_dir, "upload_Assembly")
    dest_dir = sys.argv[1]
    source_dir = sys.argv[2]
    main(dest_dir, source_dir)