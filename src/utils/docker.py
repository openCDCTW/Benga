import os
from concurrent.futures import ProcessPoolExecutor

from src.utils import operations, files


def docker_prokka_cmd(newname, outpath, inpath):
    name, ext = newname.split(".")

    args = list()
    args.append(("--cpus", "2"))
    args.append(("--outdir", "/data/" + name))
    args.append(("--prefix", name))
    prokka_ = operations.format_cmd("prokka", args, "/input/" + newname)

    args2 = list()
    args2.append(("--rm", ""))
    args2.append(("-v", outpath + ":/data"))
    args2.append(("-v", inpath + ":/input"))
    args2.append(("a504082002/prokka", ""))
    docker_prokka = operations.format_cmd("docker run", args2, prokka_)
    return docker_prokka


def prokka(newnames, outpath, inpath):
    cmds = [docker_prokka_cmd(n, outpath, inpath) for n in newnames]
    with ProcessPoolExecutor() as executor:
        executor.map(os.system, cmds)


def docker_roary_cmd(outpath, threads, identity):
    args = list()
    args.append(("-p", threads))
    args.append(("-i", identity))
    args.append(("-f", "/data/roary"))
    roary_ = operations.format_cmd("roary", args, "/data/GFF/*.gff")

    args2 = list()
    args2.append(("--rm", ""))
    args2.append(("-v", outpath + ":/data"))
    args2.append(("a504082002/roary", ""))
    docker_roary = operations.format_cmd("docker run", args2, "python /program/cmds.py " + roary_)
    return docker_roary


def roary(outpath, threads=4, ident_min=95):
    cmd = docker_roary_cmd(outpath, threads, ident_min)
    os.system(cmd)


def docker_fastx_cmd(locus_file, outpath):
    locus = os.path.splitext(locus_file)[0]

    args = list()
    args.append(("-i", files.joinpath("/data/locusfiles", locus_file)))
    args.append(("-o", files.joinpath("/data/locusfiles", locus + ".fa")))
    f = operations.format_cmd("fastx_collapser", args, "")

    args2 = list()
    args2.append(("--rm", ""))
    args2.append(("-v", outpath + ":/data"))
    args2.append(("a504082002/fastx-toolkit", ""))
    docker_fastx = operations.format_cmd("docker run", args2, f)
    return docker_fastx


def fastx(locus_files, outpath):
    cmds = [docker_fastx_cmd(file, outpath) for file in locus_files]
    with ProcessPoolExecutor() as executor:
        executor.map(os.system, cmds)
    for file in locus_files:
        os.remove(files.joinpath(outpath, "locusfiles", file))

