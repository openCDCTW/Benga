import hashlib
from functools import reduce
import functional
import uuid


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


def create_uuid():
    return str(uuid.uuid1())


def make_seqid(x):
    return hashlib.sha256(str(x).encode("ascii")).hexdigest()
