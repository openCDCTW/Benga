import os
import shutil


class ContigHandler:
    extensions = {".fna", ".fa", ".fasta"}

    def is_fasta(self, name):
        return os.path.splitext(name)[-1] in self.extensions

    def format(self, from_dir, to_dir):
        for filename in os.listdir(from_dir):
            if self.is_fasta(filename):
                new_filename = get_fileroot(filename) + ".fa"
                from_file = os.path.join(from_dir, filename)
                to_file = os.path.join(to_dir, new_filename)
                shutil.copy(from_file, to_file)


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


def get_fileroot(filename):
    return os.path.splitext(os.path.basename(filename))[0]
