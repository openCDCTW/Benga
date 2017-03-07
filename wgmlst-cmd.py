import argparse
from src.algorithms import pgdb, wgmlst, phylotree


def parse_args():
    arg_parser = argparse.ArgumentParser()

    arg_parser.add_argument(
        "-a", "--algorithm",
        required=True,
        choices=["make_db", "profiling"],
        help="Execute specified algorithm. (necessary)"
    )

    arg_parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output data directory. (necessary)"
    )

    arg_parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input data (for query) directory. (necessary)"
    )

    arg_parser.add_argument(
        "-d", "--database",
        help="Pan genome allele database for query. (necessary for profiling)"
    )

    arg_parser.add_argument(
        "-t", "--threads",
        type=int,
        default=2,
        help="Number of threads for computation. [Default: 2]",
        metavar="THREADS"
    )

    arg_parser.add_argument(
        "--identity",
        help="The minimum percentage identity for blastp. [Default: 90]",
        type=int,
        default=90,
        metavar='THRESHOLD'
    )

    arg_parser.add_argument(
        "--local",
        help="Use locally installed prokka and roary, instead of the docker version.",
        action='store_true',
        default=False
    )

    return arg_parser.parse_args()


def main():
    args = parse_args()

    input_dir = args.input
    output_dir = args.output
    db_dir = args.database
    threads = args.threads
    identity = args.identity
    docker = not args.local

    if args.algorithm == "make_db":
        pgdb.annotate_configs(input_dir, output_dir, use_docker=False)
        pgdb.make_database(output_dir, threads=threads, min_identity=identity, use_docker=docker)
    if args.algorithm == "profiling":
        wgmlst.profiling(output_dir, input_dir, db_dir,
                         threads=threads, aligcov_cut=0.5, identity=identity)


if __name__ == "__main__":
    main()
