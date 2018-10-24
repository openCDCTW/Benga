from __future__ import absolute_import, unicode_literals
from celery import shared_task
import os
from django.conf import settings
import shutil
from django.core.files import File
from src.algorithms import profiling, phylogeny
from src.utils import files
from profiling.serializers import ProfileSerializer, DendrogramSerializer
import subprocess


@shared_task
def do_profiling(batch_id, database, occr_level):
    input_dir = os.path.join(settings.MEDIA_ROOT, "sequences", batch_id)
    output_dir = os.path.join(settings.MEDIA_ROOT, "temp", batch_id)
    files.create_if_not_exist(output_dir)

    profiling.profiling(output_dir, input_dir, database, occr_level=occr_level, threads=2)

    profile_filename = os.path.join(output_dir, "wgmlst.tsv")

    dendro = phylogeny.Dendrogram()
    dendro.make_tree(profile_filename)
    newick_filename = os.path.join(output_dir, "dendrogram.newick")
    dendro.to_newick(newick_filename)
    pdf_filename = os.path.join(output_dir, "dendrogram.pdf")
    dendro.scipy_tree(pdf_filename)
    svg_filename = os.path.join(output_dir, "dendrogram.svg")
    dendro.scipy_tree(svg_filename)
    png_filename = os.path.join(output_dir, "dendrogram.png")
    dendro.scipy_tree(png_filename)
    emf_filename = os.path.join(output_dir, "dendrogram.emf")
    subprocess.call(['libreoffice', '--headless', '--convert-to', 'emf', '--outdir', output_dir, svg_filename])

    profile_data = {"id": batch_id, "file": File(open(profile_filename, "rb")),
                    "occurrence": occr_level, "database": database}
    serializer = ProfileSerializer(data=profile_data)
    if serializer.is_valid():
        serializer.save()
    else:
        print(serializer.errors)

    dendrogram_data = {"id": batch_id, "png_file": File(open(png_filename, "rb")),
                       "pdf_file": File(open(pdf_filename, "rb")), "svg_file": File(open(svg_filename, "rb")),
                       "emf_file": File(open(emf_filename, "rb")), "newick_file": File(open(newick_filename, "rb"))}
    serializer = DendrogramSerializer(data=dendrogram_data)
    if serializer.is_valid():
        serializer.save()
    else:
        print(serializer.errors)

    shutil.rmtree(output_dir)
