import argparse
import os.path
import datetime
from src.algorithms import pgdb, wgmlst, phylotree


def parse_args():
    arg_parser = argparse.ArgumentParser()

    arg_parser.add_argument(
        "-a", "--algorithm",
        required=True,
        choices=["make_db", "profiling", "tree"],
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
    docker = not args.local

    if args.algorithm == "make_db":
        pgdb.annotate_configs(input_dir, output_dir, use_docker=docker)
        pgdb.make_database(output_dir, threads=threads, use_docker=docker)
    if args.algorithm == "profiling":
        wgmlst.profiling(output_dir, input_dir, db_dir, threads=threads)
    if args.algorithm == "tree":
        dendro = phylotree.Dendrogram()
        dendro.make_tree(os.path.join(input_dir, "wgmlst.tsv"))
        date = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M")
        filename = date + "_tree"
        dendro.to_newick(os.path.join(output_dir, "{}.newick".format(filename)))
        dendro.render_on(os.path.join(output_dir, "{}.pdf".format(filename)))
        dendro.render_on(os.path.join(output_dir, "{}.svg".format(filename)))
        dendro.render_on(os.path.join(output_dir, "{}.png".format(filename)))


if __name__ == "__main__":
    main()
