
from src.utils import operations, files


def prokka(newname, outpath, inpath):
    name, ext = newname.split(".")

    args = list()
    args.append(("--outdir", files.joinpath(outpath, name)))
    args.append(("--prefix", name))
    return operations.format_cmd("prokka", args, files.joinpath(inpath, newname))


def roary(outpath, inpath, threads, ident_min):
    args = list()
    args.append(("-p", threads))
    args.append(("-i", ident_min))
    args.append(("-f", files.joinpath(outpath, "roary")))
    return operations.format_cmd("roary", args, files.joinpath(inpath, "*.gff"))


def fastx(outfile, infile):
    args = list()
    args.append(("-i", infile))
    args.append(("-o", outfile))
    return operations.format_cmd("fastx_collapser", args, "").replace("=", " ")
