import pandas as pd

from src.utils import seq


def filter_duplicates(blastp_out_file, query_length, subject_length, identity):
    '''
    blastp for locus with 95% identity, E-value < 1e-6,
    75% <= qlen/slen < 125%, 75% <= qlen/alen < 125%.
    '''
    blastp_out = pd.read_csv(blastp_out_file, sep="\t", header=None, names=seq.BLAST_COLUMNS)
    blastp_out = blastp_out[blastp_out["pident"] >= identity]
    blastp_out = blastp_out[blastp_out["qseqid"] != blastp_out["sseqid"]]
    blastp_out["qlen"] = list(map(lambda x: query_length[x], blastp_out["qseqid"]))
    blastp_out["slen"] = list(map(lambda x: subject_length[x], blastp_out["sseqid"]))
    blastp_out["qlen/slen"] = blastp_out["qlen"] / blastp_out["slen"]
    blastp_out["qlen/alen"] = blastp_out["qlen"] / blastp_out["length"]
    blastp_out = blastp_out[(0.75 <= blastp_out["qlen/slen"]) & (blastp_out["qlen/slen"] < 1.25)]
    blastp_out = blastp_out[(0.75 <= blastp_out["qlen/alen"]) & (blastp_out["qlen/alen"] < 1.25)]
    return blastp_out
