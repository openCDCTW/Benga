import os
import shutil
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


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


def replace_id(record, newid):
    return SeqRecord(Seq(str(record.seq), generic_dna), id=newid, description="")


def annotate(newname, path):
    name, ext = newname.split(".")
    cmd1 = "docker run --rm -v {path}:/data a504082002/prokka".format(**locals())
    cmd2 = "prokka --outdir /data/Assembly_ann/{name} --prefix {name} /data/Assembly/{newname}".format(**locals())
    os.system(cmd1 + " " + cmd2)


def run_fastx(locus_file, path):
    print(locus_file)
    locus, _ = locus_file.split(".")
    cmd1 = "docker run --rm -v {path}:/data sthysel/fastxtk".format(**locals())
    cmd2 = "fastx_collapser -i /data/{locus_file} -o /data/{locus}".format(**locals())
    os.system(cmd1 + " " + cmd2)


def sort_dict(d):
    a = [(k, v) for k, v in d.items()]
    return sorted(a, key=lambda x: x[0])

