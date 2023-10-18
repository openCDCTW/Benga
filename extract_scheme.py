#!/usr/bin/env python3

import os
import pandas as pd
from profiling import sequence_encoder
from utils import SQLiteDatabase


def extract(indir, outfile, threshold, locus_tag):
    database_path = outfile if os.path.splitext(outfile)[1] == '.db' else outfile + '.db'
    if os.path.exists(database_path):
        os.remove(database_path)
    pan_genome_info = os.path.join(indir, 'pan_genome_info.txt')
    pan_genome = pd.read_csv(pan_genome_info, sep='\t')
    pan_genome["locus_tag"] = [f"{locus_tag}{idx:05}" for idx, _ in enumerate(pan_genome.index, 1)]
    pan_genome.to_csv(pan_genome_info, sep='\t', index=False)

    core_genome = pan_genome[pan_genome['occurrence'] >= threshold]
    scheme = zip(core_genome.locus_tag, core_genome.dna_seq)
    alleles = (
        (sequence_encoder(dna_seq), locus_tag) 
        for locus_tag, dna_seq in zip(core_genome.locus_tag, core_genome.dna_seq)
    )
    sqlite_db = SQLiteDatabase(database_path)
    sqlite_db.create_table()
    sqlite_db.insert('scheme', scheme)
    sqlite_db.insert('alleles', alleles)
