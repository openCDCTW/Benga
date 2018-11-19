from __future__ import absolute_import, unicode_literals

import os
import pandas as pd
import shutil
import subprocess
from celery import shared_task
from django.conf import settings
from django.core.files import File

from dendrogram.serializers import DendrogramSerializer
from src.algorithms import phylogeny
from src.utils import files


@shared_task
def plot_dendrogram(batch_id):
    input_dir = os.path.join(settings.MEDIA_ROOT, "uploaded_profiles", batch_id)
    output_dir = os.path.join(settings.MEDIA_ROOT, "temp", batch_id)
    files.create_if_not_exist(output_dir)

    profiles = []
    for filename in os.listdir(input_dir):
        p = pd.read_csv(os.path.join(input_dir, filename), sep="\t", index_col=0)
        profiles.append(p)
    profiles = pd.concat(profiles, axis=1, join='inner', sort=False)

    dendro = phylogeny.Dendrogram()
    dendro.make_tree(profiles)
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

    dendrogram_data = {"id": batch_id, "png_file": File(open(png_filename, "rb")),
                       "pdf_file": File(open(pdf_filename, "rb")), "svg_file": File(open(svg_filename, "rb")),
                       "emf_file": File(open(emf_filename, "rb")), "newick_file": File(open(newick_filename, "rb"))}
    serializer = DendrogramSerializer(data=dendrogram_data)
    if serializer.is_valid():
        serializer.save()
    else:
        print(serializer.errors)

    shutil.rmtree(output_dir)
