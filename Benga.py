#!/usr/bin/env python3
import click
from profiling import profiling
from create_scheme import create
from extract_scheme import extract
from cluster import cluster


@click.group(context_settings={'help_option_names': ['-h', '--help'], 'show_default': True})
def main():
    """Bacterial Epidemiology NGs Analysis (BENGA) framework and pipeline."""


@main.command("cgmlst_profiling")
@click.option('-i', '--infile', type=click.Path(dir_okay=False), required=True, help='Path of genome.(FASTA)')
@click.option('-o', '--outfile', required=True, help='Path of output file.(TSV)')
@click.option('-d', '--database', required=True, help='Path of core-genome MLST database.')
@click.option('--prodigaltf', default="", help='Prodigal training file')
@click.option("-t", "--threads", type=int, default=1, help="Number of threads")
def cgmlst_profiling(infile, outfile, database, prodigaltf, threads):
    """Convert genomce sequence to cgMLST profile."""
    profiling(
        infile=infile, outfile=outfile, database=database, training_file=prodigaltf, threads=threads
    )


@main.command("create_scheme")
@click.option('-i', '--input-dir', type=click.Path(exists=True), required=True,
              help='Directory containing genomes.')
@click.option("-o", "--output-dir", type=click.Path(), required=True, default="benga_output",
              help="Output Directory")
@click.option("--prodigaltf", type=click.File(), default='', help='Prodigal training file')
@click.option("-t", "--threads", type=int, default=2, help="Number of threads")
def create_scheme(input_dir, output_dir, prodigaltf, threads):
    """Create wgMLST scheme"""
    create(input_dir, output_dir, threads, prodigaltf)


@main.command("extract_scheme")
@click.option('-i', '--input-dir', required=True, help='Output directory of "create" command.')
@click.option("-o", "--out-file", required=True, default='cgMLST_scheme', help="Output file of cgMLST scheme")
@click.option("-l", "--threshold", default=95, type=click.FloatRange(0, 100), help="Locus minimum occurrence")
@click.option("--locus_tag", default="Locus", help="Locus tag prefix")
def extract_scheme(input_dir, out_file, threshold, locus_tag):
    """Extract cgMLST scheme from wgMLST scheme."""
    extract(input_dir, out_file, threshold, locus_tag)


@main.command("cgmlst_cluster")
@click.option('-i', '--input-dir', required=True, help='Directory containing cgMSLT profiles.')
@click.option("-o", "--output-dir", required=True, default="cluster", help="Output Directory")
def cgmlst_cluster(input_dir, output_dir):
    """Cluster cgMLST profiles with single linkage algorithm."""
    cluster(input_dir=input_dir, output_dir=output_dir)


if __name__ == '__main__':
    main()
