import sys, hashlib, os

BUF_SIZE = 65536

ori_dir, new_dir = sys.argv[1], sys.argv[2]


def read_lines(file):
    return open(file, "r").readlines()


def compare_lines(file1, file2):
    for a, b in zip(read_lines(file1), read_lines(file2)):
        if a != b:
            print(a + " in " + file1)
            print(b + " in " + file2)
            return False
    return True


def check_sha1sum(file):
    sha1_digest = hashlib.sha1()
    with open(file, "rb") as f:
        buf = f.read(BUF_SIZE)
        while len(buf) > 0:
            sha1_digest.update(buf)
            buf = f.read(BUF_SIZE)
    checksum = sha1_digest.hexdigest()
    return checksum


def compare_sha1sum(file1, file2):
    if check_sha1sum(file1) == check_sha1sum(file2):
        return True
    else:
        return False


def check_folder(path, path2, debug=False):
    results = []
    fails = []
    compare_method = compare_lines if debug else compare_sha1sum
    for file in os.listdir(path):
        if os.path.isdir(os.path.join(path, file)):
            results.append(check_folder(os.path.join(path, file), os.path.join(path2, file)))
        else:
            result = compare_method(os.path.join(path, file), os.path.join(path2, file))
            if not result:
                fails.append(file)
            results.append(result)

    if all(results):
        print("[OK] check " + path)
        return True
    else:
        print("[!!!FAIL!!!] check " + path)
        print("The following fails:", "\t".join(fails))
        return False

check_folder(ori_dir, new_dir, True)
