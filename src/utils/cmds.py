import os.path
from src.utils import files


def form_prokka_cmd(newname, inpath, outpath):
    name, ext = newname.split(".")
    return ["prokka", "--prefix", name, "--cpus", "2",
            "--outdir", files.joinpath(outpath, name),
            files.joinpath(inpath, newname)]


def form_roary_cmd(inpath, outpath, ident_min, threads):
    return ["roary", "-p", threads, "-i", ident_min, "-s",
            "-f", files.joinpath(outpath, "roary"),
            files.joinpath(inpath, "*.gff")]


def form_prodigal_cmd(infile, outpath):
    filename = files.fasta_filename(infile)
    return ["prodigal", "-c", "-m", "-q", "-g", "11", "-i", infile,
            "-d", os.path.join(outpath, filename + ".locus.fna"),
            "-o", os.path.join(outpath, filename + ".gbk")]
