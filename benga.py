import argparse
import datetime
import os.path
import subprocess

from benga.src.algorithms import profiling, phylogeny, statistics
from src.algorithms import databases


def parse_args():
    arg_parser = argparse.ArgumentParser()

    arg_parser.add_argument(
        "-a", "--algorithm",
        required=True,
        choices=["make_db", "profiling", "tree", "statistics"],
        help="Execute specified algorithm. (necessary)"
    )

    arg_parser.add_argument(
        "-o", "--output",
        help="Output data directory. (necessary)"
    )

    arg_parser.add_argument(
        "-i", "--input",
        help="Input data (for query) directory. (necessary)"
    )

    arg_parser.add_argument(
        "-d", "--database",
        help="Pan genome allele database for query. (necessary for profiling)"
    )

    arg_parser.add_argument(
        "--occr",
        type=int,
        default=95,
        help="Level of occurrence for scheme selection. [Default: 95]",
        metavar="OCCR"
    )

    arg_parser.add_argument(
        "--drop_by_occur",
        type=float,
        default=0.0,
        help="Level of occurrence to drop. (necessary for make_db)",
        metavar="DROP"
    )

    arg_parser.add_argument(
        "-t", "--threads",
        type=int,
        default=2,
        help="Number of threads for computation. [Default: 2]",
        metavar="THREADS"
    )

    arg_parser.add_argument(
        "--no_new_alleles",
        help="Disable adding new alleles to database. [Default: False]",
        action='store_true',
        default=False
    )

    arg_parser.add_argument(
        "--no_profiles",
        help="Not generating profiles (wgmlst.tsv). This option is simply for "
             "updating database with alleles. [Default: False]",
        action='store_true',
        default=False
    )

    arg_parser.add_argument(
        "--debug",
        help="Debug mode.",
        action='store_true',
        default=False
    )

    return arg_parser.parse_args()


def main():
    args = parse_args()

    input_dir = args.input
    output_dir = args.output
    database = args.database
    occr_level = args.occr
    drop_by_occur = args.drop_by_occur
    threads = args.threads
    debug = args.debug
    new_alleles = not args.no_new_alleles
    generate_profiles = not args.no_profiles

    if args.algorithm == "make_db":
        databases.annotate_configs(input_dir, output_dir, threads=threads)
        database = databases.make_database(output_dir, drop_by_occur, threads=threads)
        statistics.calculate_loci_coverage(output_dir, output_dir, database=database)
        statistics.calculate_allele_length(output_dir, database=database)
    if args.algorithm == "statistics":
        statistics.calculate_loci_coverage(input_dir, output_dir, database=database)
        statistics.calculate_allele_length(output_dir, database=database)
    if args.algorithm == "profiling":
        profiling.profiling(output_dir, input_dir, database, threads=threads, occr_level=occr_level,
                            enable_adding_new_alleles=new_alleles, generate_profiles=generate_profiles,
                            debug=debug)
    if args.algorithm == "tree":
        dendro = phylogeny.Dendrogram()
        dendro.make_tree(os.path.join(input_dir, "wgmlst.tsv"))
        date = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M")
        filename = date + "_tree"
        dendro.to_newick(os.path.join(output_dir, "{}.newick".format(filename)))
        dendro.scipy_tree(os.path.join(output_dir, "{}.pdf".format(filename)))
        dendro.scipy_tree(os.path.join(output_dir, "{}.svg".format(filename)))
        dendro.scipy_tree(os.path.join(output_dir, "{}.png".format(filename)))
        subprocess.call(['libreoffice', '--headless', '--convert-to', 'emf',  '--outdir', output_dir,
                         os.path.join(output_dir, "{}.svg".format(filename))])


if __name__ == "__main__":
    main()
