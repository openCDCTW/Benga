import multiprocessing
import os

from src.utils import operations


class DockerFactory:
    def __init__(self):
        self._mount = dict()

    def add_mount_point(self, from_, to):
        self._mount[from_] = to

    def _key_map(self, path):
        for k, v in self._mount.items():
            if k in path:
                return path.replace(k, v)

    def map(self, paths):
        return [self._key_map(p) for p in paths]


def map_path(path, from_, to):
    return path.replace(from_, to)


def prokka_cmd(newname, outpath, inpath):
    name, ext = newname.split(".")

    args = list()
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
    pool = multiprocessing.Pool(5)
    cmds = [prokka_cmd(n, outpath, inpath) for n in newnames]
    pool.map(os.system, cmds)
    pool.close()
    pool.join()


def roary_cmd(outpath, threads, ident_min):
    args = list()
    args.append(("-p", threads))
    args.append(("-i", ident_min))
    args.append(("-f", "/data/roary"))
    roary_ = operations.format_cmd("roary", args, "/data/GFF/*.gff")

    args2 = list()
    args2.append(("--rm", ""))
    args2.append(("-v", outpath + ":/data"))
    args2.append(("a504082002/roary", ""))
    docker_roary = operations.format_cmd("docker run", args2, "python /program/cmd.py " + roary_)
    return docker_roary


def roary(outpath, threads=4, ident_min=95):
    cmd = roary_cmd(outpath, threads, ident_min)
    os.system(cmd)


def fastx_cmd(locus_file, outpath):
    locus = os.path.splitext(locus_file)[0]

    args = list()
    args.append(("-i", "/data/locusfiles/" + locus_file))
    args.append(("-o", "/data/locusfiles/" + locus + ".fa"))
    f = operations.format_cmd("fastx_collapser", args, "/data/GFF/*.gff")

    args2 = list()
    args2.append(("--rm", ""))
    args2.append(("-v", outpath + ":/data"))
    args2.append(("a504082002/fastx-toolkit", ""))
    docker_fastx = operations.format_cmd("docker run", args2, f)
    return docker_fastx


def fastx(locus_files, outpath):
    # for file in locus_files:
    #     os.system(fastx_cmd(file, outpath))
    #     print("processed " + file)
    #     os.remove("{}/locusfiles/{}".format(outpath, file))
    pool = multiprocessing.Pool(5)
    cmds = [fastx_cmd(file, outpath) for file in locus_files]
    pool.map(os.system, cmds)
    pool.close()
    pool.join()

