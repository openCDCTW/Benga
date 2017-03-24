import sys

from src.algorithms import pgdb
from src.utils import transport


def make_db(source_dir, job_dir):
    pgdb.annotate_configs(source_dir, job_dir)
    pgdb.make_database(job_dir)


def main(path):
    jobtype = "pgdb"

    receiver = transport.Receiver("127.0.0.1", 5556)
    receiver.run(path, jobtype, make_db)

if __name__ == "__main__":
    p = sys.argv[1]
    main(p)
