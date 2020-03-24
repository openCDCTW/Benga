import os.path
import subprocess

from src.utils import files


def form_prokka_cmd(genome_file, outdir, training_file):
    prefix = files.get_fileroot(genome_file)
    if os.path.exists(training_file):
        cmd = 'prokka --prefix {} --cpus 2 --outdir {} --prodigaltf {} --quiet {}'.format(prefix, outdir, training_file,
                                                                                          genome_file)
    else:
        cmd = 'prokka --prefix {} --cpus 2 --outdir {} --quiet {}'.format(prefix, outdir, genome_file)
    return cmd


def form_roary_cmd(inpath, outpath, ident_min, threads):
    cmd = "roary -p {} -i {} -s -f {} {}".format(
        threads, ident_min, os.path.join(outpath, "roary"), os.path.join(inpath, "*.gff"))
    return cmd


def form_prodigal_cmd(infile, outfile, training_file):
    if os.path.exists(training_file):
        cmd = "prodigal -c -m -q -g 11 -i {} -d {} -t {}".format(infile, outfile, training_file)
    else:
        cmd = "prodigal -c -m -q -g 11 -i {} -d {}".format(infile, outfile)
    return cmd


def execute_cmd(cmd):
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    process.communicate()
