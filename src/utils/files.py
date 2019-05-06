import os
from Bio import SeqIO
from src.utils import seq


class ContigHandler:
    def __init__(self):
        self.extensions = [".fna", ".fa", ".fasta"]
        self.seqid_template = "Contig_{}"

    def newseqid(self, j):
        return self.seqid_template.format(j)

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

    def __write_new_format(self, source_file, sink_file):
        records = []
        for j, contig in enumerate(SeqIO.parse(source_file, "fasta"), 1):
            seqid = self.newseqid(j)
            records.append(seq.new_record(seqid, str(contig.seq)))
        SeqIO.write(records, sink_file, "fasta")

    def new_format(self, from_dir, to_dir):
        filenames = []
        for filename in sorted(os.listdir(from_dir)):
            newname = self.replace_ext(filename) + ".fa"
            source_file = os.path.join(from_dir, filename)
            sink_file = os.path.join(to_dir, newname)
            self.__write_new_format(source_file, sink_file)
            filenames.append(newname)
        return filenames


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


def fasta_filename(filename):
    return os.path.basename(filename).replace(".fa", "").replace(".fna", "")

