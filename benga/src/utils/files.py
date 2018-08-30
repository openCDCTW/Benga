import os
import shutil
from Bio import SeqIO
from benga.src.utils import seq


class ContigHandler:
    def __init__(self):
        self.__namemap = {}
        self.extensions = [".fna", ".fa", ".fasta"]
        self.filename_template = "Genome_{}.fa"
        self.seqid_template = "Genome_{}::Contig_{}"

    @property
    def namemap(self):
        return self.__namemap

    def newname(self, i):
        return self.filename_template.format(i)

    def newseqid(self, i, j):
        return self.seqid_template.format(i, j)

    def isfasta(self, name):
        for ext in self.extensions:
            if name.endswith(ext):
                return True
        else:
            return False

    def replace_ext(self, xs):
        ys = xs
        for ext in self.extensions:
            ys = ys.replace(ext, "")
        return ys

    def __write_new_format(self, source_file, sink_file, i):
        records = []
        for j, contig in enumerate(SeqIO.parse(source_file, "fasta"), 1):
            seqid = self.newseqid(i, j)
            records.append(seq.new_record(seqid, str(contig.seq)))
        SeqIO.write(records, sink_file, "fasta")

    def new_format(self, from_dir, to_dir, replace_ext=True):
        for i, filename in enumerate(sorted(os.listdir(from_dir)), 1):
            newname = self.newname(i)
            if replace_ext:
                self.__namemap[self.replace_ext(newname)] = self.replace_ext(filename)
            else:
                self.__namemap[newname] = filename
            source_file = os.path.join(from_dir, filename)
            sink_file = os.path.join(to_dir, newname)
            self.__write_new_format(source_file, sink_file, i)


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
        os.makedirs(path)


def clear_folder(path):
    if os.path.exists(path):
        shutil.rmtree(path)
    os.mkdir(path)


def recursive_chown(dir, user):
    for root, dirs, files in os.walk(dir):
        for d in dirs:
            shutil.chown(joinpath(root, d), user)
        for f in files:
            shutil.chown(joinpath(root, f), user)


def fasta_filename(filename):
    return os.path.basename(filename).replace(".fa", "").replace(".fna", "")

