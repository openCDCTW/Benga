import os
import shutil


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


def recursive_chown(dir, user):
    for root, dirs, files in os.walk(dir):
        for d in dirs:
            shutil.chown(joinpath(root, d), user)
        for f in files:
            shutil.chown(joinpath(root, f), user)


def replace_ext(x, extensions=[".fna", ".fa"]):
    y = x
    for ext in extensions:
        y = y.replace(ext, "")
    return y


def fasta_filename(filename):
    return os.path.basename(filename).replace(".fa", "").replace(".fna", "")

