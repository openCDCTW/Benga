import datetime
import os.path
import json
import psycopg2
from src.utils import db
from src.algorithms import wgmlst, phylotree

INDIR = "input"
OUTDIR = "output"


def profiling_api(args):
    batch_id, database, occr_level = args
    input_dir = os.path.join(INDIR, batch_id)
    output_dir = os.path.join(OUTDIR, batch_id)
    wgmlst.profiling(output_dir, input_dir, database, occr_level=occr_level, threads=2)
    profile_created = datetime.datetime.now()

    with open(os.path.join(output_dir, "namemap.json"), "r") as file:
        names = json.loads(file.read())
    profile_filename = os.path.join(output_dir,
                                    "cgMLST_{}_{}_{}.tsv".format(database, occr_level, batch_id[0:8]))
    dendro = phylotree.Dendrogram()
    dendro.make_tree(profile_filename, names)
    dendro_created = datetime.datetime.now()
    newick_filename = os.path.join(output_dir, "dendrogram_{}.newick".format(batch_id[0:8]))
    dendro.to_newick(newick_filename)
    pdf_filename = os.path.join(output_dir, "dendrogram_{}.pdf".format(batch_id[0:8]))
    dendro.scipy_tree(pdf_filename)
    svg_filename = os.path.join(output_dir, "dendrogram_{}.svg".format(batch_id[0:8]))
    dendro.scipy_tree(svg_filename)
    png_filename = os.path.join(output_dir, "dendrogram_{}.png".format(batch_id[0:8]))
    dendro.scipy_tree(png_filename)

    with open(profile_filename, "rb") as file:
        profile_file = file.read()
    sql = "INSERT INTO profile (id,created,file,occurrence,database) VALUES(%s,%s,%s,%s,%s);"
    data = (batch_id, profile_created, psycopg2.Binary(profile_file), occr_level, database)
    db.to_sql(sql, data, database="profiling")

    with open(png_filename, "rb") as file:
        png_file = file.read()
    with open(pdf_filename, "rb") as file:
        pdf_file = file.read()
    with open(svg_filename, "rb") as file:
        svg_file = file.read()
    with open(newick_filename, "rb") as file:
        newick_file = file.read()
    sql = "INSERT INTO dendrogram (id,created,png_file,pdf_file,svg_file,newick_file) VALUES(%s,%s,%s,%s,%s);"
    data = (batch_id, dendro_created, psycopg2.Binary(png_file.read()),
            psycopg2.Binary(pdf_file.read()), psycopg2.Binary(svg_file.read()),
            psycopg2.Binary(newick_file.read()))
    db.to_sql(sql, data, database="profiling")
