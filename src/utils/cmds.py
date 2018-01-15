import os.path
import subprocess
from src.utils import files


def form_prokka_cmd(newname, inpath, outpath):
    name, ext = newname.split(".")
    args = ["prokka", "--prefix", name, "--cpus", "2", "--outdir", files.joinpath(outpath, name),
            files.joinpath(inpath, newname)]
    return " ".join(args)


def form_roary_cmd(inpath, outpath, ident_min, threads):
    args = ["roary", "-p", threads, "-i", ident_min, "-s", "-f", files.joinpath(outpath, "roary"),
            files.joinpath(inpath, "*.gff")]
    return " ".join(args)


def form_prodigal_cmd(infile, outpath):
    filename = files.fasta_filename(infile)
    args = ["prodigal", "-c", "-m", "-q", "-g", "11", "-i", infile,
            "-d", os.path.join(outpath, filename + ".locus.fna"),
            "-o", os.path.join(outpath, filename + ".gbk")]
    return " ".join(args)


def execute_cmd(args):
    cmd = args
    subprocess.run(cmd, shell=True)