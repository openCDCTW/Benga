import json
from collections import defaultdict, OrderedDict
import functional

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import fastcluster
from scipy.cluster.hierarchy import dendrogram, to_tree
from Bio import SeqIO

from utils import *


namemap = None
SCHEMES = ["locusAP_core", "locusAP_core_10p_dispensable", "locusAP_core_20p_dispensable",
           "locusAP_core_30p_dispensable", "locusAP_core_40p_dispensable", "locusAP_core_50p_dispensable",
           "locusAP_core_60p_dispensable", "locusAP_core_70p_dispensable", "locusAP_core_80p_dispensable",
           "locusAP_core_90p_dispensable", "locusAP_core_100p_dispensable", "locusAP_pan",
           "locusAP_unique", "wgMLST_core", "wgMLST_core_10p_dispensable",
           "wgMLST_core_20p_dispensable", "wgMLST_core_30p_dispensable", "wgMLST_core_40p_dispensable",
           "wgMLST_core_50p_dispensable", "wgMLST_core_60p_dispensable", "wgMLST_core_70p_dispensable",
           "wgMLST_core_80p_dispensable", "wgMLST_core_90p_dispensable", "wgMLST_core_100p_dispensable",
           "wgMLST_pan", "wgMLST_unique"]


def annotate_configs(temp_dir, source_dir, output_dir, number_per_docker=10):
    print("## 02_reformatContigFile ##")

    filenames = parse_filenames(source_dir)
    working_dir = joinpath(temp_dir, "Assembly")
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

    annotate_dir = joinpath(temp_dir, "Assembly_ann")
    create_if_not_exist(annotate_dir)
    prokka(namemap.values())

    print("--Finished")


    print("## 04_make_ffn-gff ##")

    gff_dir = joinpath(temp_dir, "GFF")
    ffn_dir = joinpath(temp_dir, "FFN")
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

    with open(joinpath(output_dir, "nonCDS.json"), "w") as file:
        file.write(json.dumps(noncds))

    print("--Finished")


def make_profiles(temp_dir, output_dir, threads=2, ident_min=95, number_per_docker=500):
    print("## 01_reformat_ffn ##")

    ffn_dir = joinpath(temp_dir, "FFN")
    ffn2_dir = joinpath(temp_dir, "FFN2")
    create_if_not_exist(ffn2_dir)

    for filename in os.listdir(ffn_dir):
        handle = SeqIO.parse(joinpath(ffn_dir, filename), "fasta")
        records = [new_record(record.id, str(record.seq)) for record in handle]
        SeqIO.write(records, joinpath(ffn2_dir, filename), "fasta")

    print("--Finished")


    # calculate the pan genome
    print("## 02_runRoary ##")

    cmd = "roary -e --mafft -p {threads} -i {ident_min} -f /data/roary /data/GFF/*.gff".format(**locals())
    os.system(cmd)

    print("--Finished")


    print("## 03_parse_csvfile ##")

    roary_dir = joinpath(temp_dir, "roary")
    locus_dir = joinpath(temp_dir, "locusfiles")
    create_if_not_exist(locus_dir)

    roary_result = pd.read_csv(joinpath(roary_dir, "gene_presence_absence.csv"))
    roary_result.index = pd.Index(map(lambda x: "SAL{0:07d}".format(x+1), roary_result.index), name="locus")

    mapping = roary_result[roary_result["No. isolates"] == roary_result["No. sequences"]]
    mapping = mapping[["Gene", "No. isolates", "No. sequences", "Annotation"]]
    mapping.to_csv(joinpath(temp_dir, "locusmapping.txt"), sep="\t")

    paralog = roary_result[roary_result["No. isolates"] != roary_result["No. sequences"]]
    paralog = paralog[["Gene", "No. isolates", "No. sequences", "Annotation"]]
    paralog.to_csv(joinpath(temp_dir, "paralog.txt"), sep="\t", index=False)

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
    fastx(locus_files)

    print("--Finished")


    print("## 05_renumberAlleles ##")

    for locus_file in locus_files:
        file_path = joinpath(locus_dir, locus_file.split(".")[0])

        records = functional.seq(
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
        name = joinpath(temp_dir, "locusfiles", filename)
        seq_len = [(str(record.seq), len(str(record.seq))) for record in SeqIO.parse(name, "fasta")]
        maxlen_seq, maxlen = max(seq_len, key=lambda x: x[1])
        refseqs[filename] = maxlen_seq

    records = [new_record(key, value) for key, value in refseqs.items()]
    SeqIO.write(records, joinpath(output_dir, "panRefSeq.fa"), "fasta")

    with open(joinpath(output_dir, "panRefSeq.json"), "w") as file:
        file.write(json.dumps(refseqs))

    print("--Finished")


    print("## 07_makeSchems ##")

    scheme_dir = joinpath(temp_dir, "scheme")
    create_if_not_exist(scheme_dir)

    total_isolates = len(os.listdir(ffn2_dir))
    dispatcher = defaultdict(list)
    for index, row in pd.read_csv(joinpath(temp_dir, "locusmapping.txt"), sep="\t").iterrows():
        locus = row["locus"]
        occ = row["No. isolates"]
        occrate = occ / total_isolates

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


    print("## 08_output ##")

    shutil.copy(joinpath(temp_dir, "roary", "summary_statistics.txt"), output_dir)
    shutil.copy(joinpath(temp_dir, "locusmapping.txt"), output_dir)
    shutil.move(scheme_dir, output_dir)
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


def query_blast(m, dest_dir, dbsource, query, aligcov_cut, pident_cut, seqlen, threads, cols):
    blast_out = joinpath(dest_dir, "blast.out")
    locidb = joinpath(dest_dir, "lociDB")

    os.system(makeblastdb(dbsource, "nucl", locidb))
    os.system(blastn(threads, query, locidb, blast_out, "'6 {items}'".format(items=" ".join(cols))))

    result = pd.read_csv(blast_out, sep="\t", header=None, names=cols)
    result["qlen"] = list(map(lambda x: seqlen[x], result["qseqid"]))
    result["aligcov"] = (result["length"] - result["gapopen"]) / result["qlen"]
    result = result[result["aligcov"] >= aligcov_cut][result["pident"] >= pident_cut]
    for index, row in result.iterrows():
        row_elements = [row["aligcov"], row["length"], row["qlen"], row["pident"], row["sseqid"], row["sstart"],
                        row["send"]]
        m[row["qseqid"]] = "\t".join(map(str, row_elements))
    os.remove(blast_out)
    return m


def maxlen_locus(locus, pair):
    if all(map(lambda x: x == 0, pair)):
        return locus, 0
    else:
        # only access the last of max length allele number in values
        return locus, max(reversed(pair), key=lambda x: x[1])[0]


def wg_profiling(temp_dir, input_dir, blastn_threads=2, aligcov_cut=0.5, pident_cut=90):
    print("## 01_reformatContigFile ##")

    assemble_dir = joinpath(temp_dir, "query_assembly")
    create_if_not_exist(assemble_dir)

    global namemap
    namemap = {}
    for i, filename in enumerate(sorted(os.listdir(input_dir)), 1):
        file = SeqIO.parse(joinpath(input_dir, filename), "fasta")
        records = []
        for j, record in enumerate(file, 1):
            newid = "Assembly_{i}::Contig_{j}".format(**locals())
            records.append(new_record(newid, str(record.seq)))

        newname = "Assembly_{i}.fa".format(**locals())
        SeqIO.write(records, joinpath(assemble_dir, newname), "fasta")
        namemap[newname] = filename

    print("--Finished")

    print("## 02_locusProfiling ##")

    scheme_dir = joinpath(temp_dir, "scheme")
    blastquery = joinpath(scheme_dir, "pan.txt")
    shutil.copytree(joinpath(input_dir, "scheme"), scheme_dir)

    seqlen = {record.id: len(record.seq) for record in SeqIO.parse(blastquery, "fasta")}
    m = [(id, "") for id in seqlen.keys()]
    m = OrderedDict(sorted(m, key=lambda x: x[0]))

    blast_cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend",
                  "sstart", "send", "evalue", "bitscore"]
    for filename in os.listdir(assemble_dir):
        m_ = query_blast(m.copy(), temp_dir, joinpath(assemble_dir, filename), blastquery,
                         aligcov_cut, pident_cut, seqlen, blastn_threads, blast_cols)
        output = functional.pseq(
            m_.items()
        ).map(
            lambda k, v: (k, 1 if len(v) > 1 else 0)
        ).map(
            lambda x, hit: x + "\t" + hit
        ).to_list()

        name = filename.split(".")[0]
        with open(joinpath(temp_dir, "locusAP." + name + ".list"), "w") as file:
            file.write("\n".join(output))

    print("--Finished")

    print("## 03_alleleAllocate ##")

    for filename in os.listdir(assemble_dir):
        qfn2 = filename.split(".")[0]
        print(qfn2)
        os.mkdir(joinpath(temp_dir, qfn2))

        blastinfile = joinpath(assemble_dir, filename)
        os.system(makeblastdb(blastinfile, "nucl", joinpath(temp_dir, "AssemblyDB_" + qfn2)))

        allprofile = defaultdict(list)
        newloc = []
        aldic = {}
        for line in open(joinpath(temp_dir, "locusAP." + qfn2 + ".list"), "r").read().splitlines():
            locus, hits = line.split("\t")
            if hits == "0":
                allprofile[locus] = [0]
                newloc.append(locus)
            else:
                qfile = joinpath(input_dir, "locusfiles", locus + ".new")
                for record in SeqIO.parse(qfile, "fasta"):
                    alleleno = record.id.split("::")[1]
                    aldic[(locus, alleleno)] = len(record.seq)
                os.system("cat {qfile} >> {temp_dir}/qfile.fa".format(**locals()))

        os.system(blastn(blastn_threads, joinpath(temp_dir, "qfile.fa"), joinpath(temp_dir, "AssemblyDB_" + qfn2),
                         joinpath(temp_dir, "blast." + qfn2 + ".out"), "'6 {blast_items}'".format(**locals())))
        os.remove(joinpath(temp_dir, "qfile.fa"))
        blast = pd.read_csv(joinpath(temp_dir, "blast." + qfn2 + ".out"), sep="\t", header=None, names=blast_cols)
        blast["loc"] = blast["qseqid"].str.split("::").str[0]
        blast["allno"] = blast["qseqid"].str.split("::").str[1]

        for _, row in blast[(blast["pident"] == 100) & (blast["mismatch"] == 0) & (blast["gapopen"] == 0)].iterrows():
            if row["length"] == aldic[(row["loc"], row["allno"])]:
                allprofile[row["loc"]].append((row["allno"], row["length"]))

        # New alleles
        new_alleles = []
        records = []
        assdic = {record.id: str(record.seq) for record in SeqIO.parse(blastinfile, "fasta")}
        for _, row in blast[~blast["loc"].isin(allprofile.keys())].iterrows():
            allprofile[row["loc"]] = [0]
            new_alleles.append(row["loc"])
            records.append(new_record(row["loc"], assdic[row["sseqid"]]))

        SeqIO.write(drop_duplicate(records, lambda x: x.id), joinpath(temp_dir, qfn2, "tot.1.new"), "fasta")
        with open(joinpath(temp_dir, qfn2, "tot.new"), "w") as file:
            file.write("\n".join(drop_duplicate(new_alleles)))

        output = functional.pseq(
            allprofile.items()
        ).map(
            lambda locus, pair: maxlen_locus(locus, pair)
        ).map(
            lambda a: sorted(a, key=lambda x: x[0])
        ).map(
            lambda locus, pair: str(locus) + "\t" + str(pair)
        ).to_list()
        with open(joinpath(temp_dir, "allele." + qfn2 + ".profile"), "w") as file:
            file.write("\n".join(output))

    print("--Finished")

    print("## 04_catchSchemeProfile ##")

    scheme = {}
    for line in open(joinpath(scheme_dir, "scheme.txt"), "r").read().splitlines():
        token = line.split("\t")
        scheme[token[0]] = token[1:]
        os.mkdir(joinpath(temp_dir, "locusAP_" + token[0]))
        os.mkdir(joinpath(temp_dir, "wgMLST_" + token[0]))

    for filename in os.listdir(assemble_dir):
        ass = filename.split(".")[0]
        for key, value in scheme.items():
            locusAP_output = functional.pseq(
                open(joinpath(temp_dir, "locusAP." + ass + ".list"), "r").read().splitlines()
            ).map(
                lambda x: x.split("\t")
            ).filter(
                lambda x: x[0] in value
            ).map(
                lambda x: x[0] + "\t" + x[1]
            ).to_list()

            wgMLST_output = functional.pseq(
                open(joinpath(temp_dir, "allele." + ass + ".profile"), "r").read().splitlines()
            ).map(
                lambda x: x.split("\t")
            ).filter(
                lambda x: x[0] in value
            ).map(
                lambda x: x[0] + "\t" + x[1]
            ).to_list()

            with open("{temp_dir}/locusAP_{key}/locusAP.{ass}.{key}.list".format(**locals()), "w") as file:
                file.write("\n".join(locusAP_output))
            with open("{temp_dir}/wgMLST_{key}/allele.{ass}.{key}.list".format(**locals()), "w") as file:
                file.write("\n".join(wgMLST_output))

    print("--Finished")

    print("## 05_outformat ##")

    locusap_dir = joinpath(scheme_dir, "locusAP")
    clear_folder(locusap_dir)

    wgmlst_dir = joinpath(scheme_dir, "wgMLST")
    clear_folder(wgmlst_dir)

    os.system("mv {temp_dir}/locusAP_* {locusap_dir}".format(**locals()))
    os.system("mv {temp_dir}/wgMLST_* {wgmlst_dir}".format(**locals()))

    dir_names = ["_core", "_pan", "_core_10p_dispensable", "_core_20p_dispensable", "_core_30p_dispensable",
                 "_core_40p_dispensable", "_core_50p_dispensable", "_core_60p_dispensable", "_core_70p_dispensable",
                 "_core_80p_dispensable", "_core_90p_dispensable", "_core_100p_dispensable"]
    for name in dir_names:
        locusap_name = joinpath(locusap_dir, "locusAP" + name)
        wgmlst_name = joinpath(wgmlst_dir, "wgMLST" + name)
        os.system("ls {locusap_name} | grep 'locusAP'> " + joinpath(locusap_name, "files.list"))
        os.system("ls {wgmlst_name} | grep 'allele'> " + joinpath(wgmlst_name, "files.list"))

    os.system("rm -rf Assembly*")

    print("--Finished")


def getNewick(node, newick, parentdist, leaf_names):
    if node.is_leaf():
        return "{}:{:.2f}{}".format(leaf_names[node.id], parentdist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):{:.2f}{}".format(parentdist - node.dist, newick)
        else:
            newick = ");"
        newick = getNewick(node.get_left(), newick, node.dist, leaf_names)
        newick = getNewick(node.get_right(), ",{}".format(newick), node.dist, leaf_names)
        newick = "({}".format(newick)
        return newick


def clustering(temp_dir, output_dir, scheme_selection):
    # makeDiffMatrix
    scheme_path = {s: joinpath("scheme", s.split("_")[0], s) for s in SCHEMES}
    target_dir = joinpath(temp_dir, scheme_path[scheme_selection])

    global namemap

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
        for line in open(joinpath(target_dir, filename)).read().splitlines():
            if line == "":
                continue
            locus, hit = line.split("\t")
            column[locus] = True if hit != 0 else False
        sample = filename.split(".")[1]
        matrix[sample] = column

    # save profiles
    matrix.to_csv(joinpath(output_dir, "profiles.tsv"), sep="\t")

    # build phylogenetic tree
    Y = fastcluster.linkage(matrix, method="average", metric="hamming")

    # save to pdf
    fig = plt.figure()
    dendrogram(Y, orientation="left")
    pp = PdfPages(joinpath(output_dir, "phylogenetic_tree.pdf"))
    pp.savefig(fig)
    plt.close()

    # save to newick
    tree = to_tree(Y, False)
    newick = getNewick(tree, "", tree.dist, leaf_names=matrix.columns)
    with open(joinpath(output_dir, "phylogenetic_tree.newick"), "w") as file:
        file.write(newick)


def make_database(temp_dir, source_dir, output_dir):
    annotate_configs(temp_dir, source_dir, output_dir)
    make_profiles(temp_dir, output_dir)


def build_tree(database_dir, query_dir, output_dir, scheme_selection):
    wg_profiling(database_dir, query_dir)
    clustering(database_dir, output_dir, scheme_selection)


if __name__ == "__main__":
    # dest_dir = sys.argv[1]
    # source_dir = sys.argv[2]
    scheme_selection = "locusAP_core_100p_dispensable"
    # make_database("/data", "/input", "/output")
    build_tree(dest_dir, scheme_selection)

