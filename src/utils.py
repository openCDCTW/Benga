import os
import shutil
import multiprocessing
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


def use_docker(docker, path, procs=2):
    def method_decorator(func):
        def wrapper(items):
            cmd = "docker run --rm -v {path}:/data {docker} python /program/batch.py {procs}".format(
                **locals())
            cmds = " ".join(map(func, items))
            os.system(cmd + " " + cmds)
        return wrapper
    return method_decorator


def use_parallel(procs=4):
    def method_decorator(func):
        def wrapper(items):
            pool = multiprocessing.Pool(procs)
            cmds = list(map(func, items))
            pool.map(os.system, cmds)
            pool.close()
            pool.join()
        return wrapper
    return method_decorator


@use_parallel(4)
def prokka(newname):
    name, ext = newname.split(".")
    return "prokka --outdir /data/Assembly_ann/{name} --prefix {name} /data/Assembly/{newname}".format(**locals())


@use_parallel(4)
def fastx(locus_file):
    locus, _ = locus_file.split(".")
    return "fastx_collapser -i /data/locusfiles/{locus_file} -o /data/locusfiles/{locus}".format(**locals())


def sort_dict(d):
    a = [(k, v) for k, v in d.items()]
    return sorted(a, key=lambda x: x[0])


def makeblastdb(input, dbtype, output):
    return "makeblastdb -in {input} -dbtype {dbtype} -out {output}".format(**locals())


def blastn(num_threads, query, db, output, outfmt):
    return "blastn -num_threads {num_threads} -query {query} -db {db} -out {output} -outfmt {outfmt}".format(**locals())


