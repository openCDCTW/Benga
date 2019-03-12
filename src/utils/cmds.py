import os.path
import subprocess

from src.utils import files

DIR_PATH = os.path.dirname(os.path.realpath(__file__))
BINARIES_PATH = os.path.abspath(os.path.join(DIR_PATH, "..", "..", "binaries", "linux"))
MODELS_PATH = os.path.abspath(os.path.join(DIR_PATH, "..", "..", 'models'))


def form_prokka_cmd(newname, inpath, outpath):
    name, ext = newname.split(".")
    args = ["prokka", "--prefix", name, "--cpus", "2", "--outdir", os.path.join(outpath, name),
            os.path.join(inpath, newname)]
    return " ".join(map(str, args))


def form_roary_cmd(inpath, outpath, ident_min, threads):
    args = ["roary", "-p", threads, "-i", ident_min, "-s", "-f", os.path.join(outpath, "roary"),
            os.path.join(inpath, "*.gff")]
    return " ".join(map(str, args))


def form_prodigal_cmd(infile, outpath, model):
    filename = files.fasta_filename(infile)
    args = [os.path.join(BINARIES_PATH, "prodigal"), "-c", "-m", "-q", "-g", "11",
            "-t", os.path.join(MODELS_PATH, model + '.trn'),
            "-i", infile,
            "-d", os.path.join(outpath, filename + ".locus.fna"),
            "-o", os.path.join(outpath, filename + ".gbk")]
    return " ".join(map(str, args))


def execute_cmd(args):
    cmd = args
    subprocess.run(cmd, shell=True)
