#!/usr/bin/env python3

import os
import argparse
import pandas as pd
from utils import SQLiteDatabase


def main():
    parser = argparse.ArgumentParser(description="Defining the cgMLST scheme")
    parser.add_argument(
        "-i", "--input", required=True, help="Path of future create_scheme.py output"
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Path of sqlite database"
    )
    parser.add_argument(
        "-l", "--threshold", default=95, type=int,
        help="Locus minimum occurrence (0-100)"
    )
    parser.add_argument(
        "--locus_tag", default="Locus",
        help="Locus tag prefix (default = Locus)"
    )
    args = parser.parse_args()
    database_path = args.output if os.path.splitext(args.output)[1] == '.db' else args.output + '.db'
    locus_tag = args.locus_tag

    if os.path.exists(database_path):
        os.remove(database_path)
    pan_genome_info = os.path.join(args.input, 'pan_genome_info.txt')
    pan_genome = pd.read_csv(pan_genome_info, sep='\t')
    pan_genome["locus_tag"] = [f"{locus_tag}{idx:05}" for idx, _ in enumerate(pan_genome.index, 1)]
    pan_genome.to_csv(pan_genome_info, sep='\t', index=False)

    core_genome = pan_genome[pan_genome['occurrence'] >= args.threshold]
    values = zip(core_genome.locus_tag, core_genome.dna_seq)
    sqlite_db = SQLiteDatabase(database_path)
    sqlite_db.create_table()
    sqlite_db.insert('scheme', values)


if __name__ == '__main__':
    main()
