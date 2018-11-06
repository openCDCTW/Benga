import click
import datetime
import os.path
import subprocess
from src.algorithms import databases, profiling, phylogeny, statistics


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.group(context_settings=CONTEXT_SETTINGS)
def main():
    """Welcome to Bacterial Epidemiology NGs Analysis (BENGA) pipeline."""


@main.command("makedb", short_help="Make pan-genome allele database",
              context_settings=CONTEXT_SETTINGS)
@click.option('--drop_by_occur', default=0.0, metavar="<float>", type=float,
              help="Level of occurrence to drop.")
@click.option('-t', '--threads', default=8, metavar="<int>", type=int,
              help="Number of threads for computation. [Default: 8]")
@click.argument('input_dir', type=click.Path(exists=True))
@click.argument('output_dir', type=click.Path(exists=True))
def makedb(input_dir, output_dir, drop_by_occur, threads):
    """Make database with fasta files in INPUT_DIR and output accessory results in OUTPUT_DIR."""
    databases.annotate_configs(input_dir, output_dir, threads=threads)
    database = databases.make_database(output_dir, drop_by_occur, threads=threads)
    statistics.calculate_loci_coverage(output_dir, output_dir, database=database)
    statistics.calculate_allele_length(output_dir, database=database)


@main.command("stats", short_help="Make database statistics",
              context_settings=CONTEXT_SETTINGS)
@click.argument('database', type=str)
@click.argument('input_dir', type=click.Path(exists=True))
@click.argument('output_dir', type=click.Path(exists=True))
def statistics(input_dir, output_dir, database):
    """Plot statistics of DATABASE and profile from INPUT_DIR, and then output to OUTPUT_DIR."""
    statistics.calculate_loci_coverage(input_dir, output_dir, database=database)
    statistics.calculate_allele_length(output_dir, database=database)


@main.command("profile", short_help="Make profiles against database",
              context_settings=CONTEXT_SETTINGS)
@click.option('-o', '--occrrence', default=95, metavar="<int>", type=int,
              help="Level of occurrence for scheme selection. [Default: 95]")
@click.option('-t', '--threads', default=8, metavar="<int>", type=int,
              help="Number of threads for computation. [Default: 8]")
@click.option('--not-extend', default=False, is_flag=True,
              help="Disable allele extension. [Default: Enable]")
@click.option('--no-profiles', default=False, is_flag=True,
              help="Disable generating profiles (profile.tsv). [Default: Enable]")
@click.option('--debug', default=False, is_flag=True,
              help="Print additional information.")
@click.argument('database', type=str)
@click.argument('input_dir', type=click.Path(exists=True))
@click.argument('output_dir', type=click.Path(exists=True))
def profiling(input_dir, output_dir, database, threads, occrrence, not_extend, no_profiles,
              debug):
    """Make profiles with fasta files in INPUT_DIR against DATABASE, and then output to OUTPUT_DIR."""
    # profiling.profiling(output_dir, input_dir, database, threads=threads, occr_level=occrrence,
    #                     enable_adding_new_alleles=(not not_extend), generate_profiles=(not no_profiles),
    #                     debug=debug)


@main.command("tree", short_help="Plot dendrogram",
              context_settings=CONTEXT_SETTINGS)
@click.argument('input_dir', type=click.Path(exists=True))
@click.argument('output_dir', type=click.Path(exists=True))
def tree(input_dir, output_dir):
    """Plot dendrogram with profile.tsv file in INPUT_DIR, and output to OUTPUT_DIR."""
    dendro = phylogeny.Dendrogram()
    dendro.make_tree(os.path.join(input_dir, "profile.tsv"))
    date = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M")
    filename = date + "_tree"
    dendro.to_newick(os.path.join(output_dir, "{}.newick".format(filename)))
    dendro.scipy_tree(os.path.join(output_dir, "{}.pdf".format(filename)))
    dendro.scipy_tree(os.path.join(output_dir, "{}.svg".format(filename)))
    dendro.scipy_tree(os.path.join(output_dir, "{}.png".format(filename)))
    subprocess.call(['libreoffice', '--headless', '--convert-to', 'emf', '--outdir', output_dir,
                     os.path.join(output_dir, "{}.svg".format(filename))])


if __name__ == "__main__":
    main()
