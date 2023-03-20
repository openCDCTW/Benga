import argparse


def main():
    parser = argparse.ArgumentParser(description="Bacterial Epidemiology NGs Analysis (BENGA) framework and pipeline.")
    subparsers = parser.add_subparsers()

    profiling_parser = subparsers.add_parser("profiling", help="Convert genomce sequence to cgMLST profile.")
    profiling_parser.add_argument(
        "-i", "--input",
        required=True,
        help="Path of query genome."
    )
    profiling_parser.add_argument(
        "-o", "--output",
        required=True,
        help="Path of output file."
    )
    profiling_parser.add_argument(
        "-d", "--database",
        required=True,
        help="Path of core-genome MLST database."
    )
    profiling_parser.add_argument(
        "--prodigaltf",
        default='',
        help="Path of prodigal training file. default: None"
    )
    profiling_parser.add_argument(
        "-t", "--threads",
        type=int,
        default=1,
        help="Number of threads. default: 1"
    )

    create_scheme_parser = subparsers.add_parser("create_scheme", help="Create wgMLST scheme.")
    create_scheme_parser.add_argument(
        "-i", "--input-files", nargs='+',
        required=True,
        help="Path of query assembly files."
    )
    create_scheme_parser.add_argument(
        "-o", "--output-dir",
        required=True,
        help="Path of output directory."
    )
    create_scheme_parser.add_argument(
        "--prodigaltf",
        default='',
        help="Path of prodigal training file. default:''"
    )
    create_scheme_parser.add_argument(
        "-t", "--threads",
        type=int,
        default=2,
        help="Number of threads. default: 2"
    )

    extract_scheme_parser = subparsers.add_parser("extract_scheme", help="Extract cgMLST scheme from wgMLST scheme.")
    extract_scheme_parser.add_argument(
        "-i", "--input", required=True, help="Path of future create_scheme.py output"
    )
    extract_scheme_parser.add_argument(
        "-o", "--output", required=True, help="Path of sqlite database"
    )
    extract_scheme_parser.add_argument(
        "-l", "--threshold", default=95, type=int,
        help="Locus minimum occurrence (0-100)"
    )
    extract_scheme_parser.add_argument(
        "--locus_tag", default="Locus",
        help="Locus tag prefix (default = Locus)"
    )
    args = parser.parse_args()


if __name__ == '__main__':
    main()
