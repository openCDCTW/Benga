import os
import shutil
from functools import reduce
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def depart(x, y, number):
    if isinstance(x, list):
        if len(x[-1]) == number:
            x.append([y])
        else:
            x[-1].append(y)
        return x
    else:
        return [[x, y]]


def partition(lst, number):
    return reduce(lambda x, y: depart(x, y, number), lst)


def joinpath(a, *args):
    return os.path.join(a, *args)


def drop_duplicate(l, idfun=None):
   if idfun is None:
       idfun = lambda x: x
   seen = {}
   result = []
   for item in l:
       marker = idfun(item)
       if marker not in seen:
           seen[marker] = 1
           result.append(item)
   return result


def create_if_not_exist(path):
    if not os.path.exists(path):
        os.mkdir(path)


def clear_folder(path):
    if os.path.exists(path):
        shutil.rmtree(path)
    os.mkdir(path)


def parse_filenames(path, ext=".fna"):
    return [name for name in os.listdir(path) if name.endswith(ext)]


def new_record(seqid, seq, desc=""):
    return SeqRecord(Seq(seq, generic_dna), id=seqid, description=desc)


def replace_id(record, newid):
    return new_record(newid, str(record.seq))


def form_prokka_cmd(newname):
    name, ext = newname.split(".")
    return "'prokka --outdir /data/Assembly_ann/{name} --prefix {name} /data/Assembly/{newname}'".format(**locals())


def annotate(newnames, path, procs=2):
    cmd = "docker run --rm -v {path}:/data a504082002/prokka python /program/batch.py {procs}".format(**locals())
    cmds = " ".join(map(form_prokka_cmd, newnames))
    os.system(cmd + " " + cmds)


def form_fastx_cmd(locus_file):
    locus, _ = locus_file.split(".")
    return "'fastx_collapser -i /data/{locus_file} -o /data/{locus}'".format(**locals())


def run_fastx(locus_files, path, procs=2):
    print(locus_files)
    cmd = "docker run --rm -v {path}:/data a504082002/fastx-toolkit python /program/batch.py {procs}".format(**locals())
    cmds = " ".join(map(form_fastx_cmd, locus_files))
    os.system(cmd + " " + cmds)


def sort_dict(d):
    a = [(k, v) for k, v in d.items()]
    return sorted(a, key=lambda x: x[0])

