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


def read_profiles(input_dir):
    files = list(filter(lambda x: x.endswith(".tsv"), os.listdir(input_dir)))
    if len(files) == 1:
        profiles = pd.read_csv(os.path.join(input_dir, files[0]), sep="\t", index_col=0)
    else:
        profiles = []
        for filename in files:
            p = pd.read_csv(os.path.join(input_dir, filename), sep="\t", index_col=0)
            profiles.append(p)
        profiles = pd.concat(profiles, axis=1, join='inner', sort=False)
    return profiles


def plot(input_dir, output_dir):
    profiles = read_profiles(input_dir)
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
    return emf_filename, newick_filename, pdf_filename, png_filename, svg_filename


def save(batch_id, emf_filename, newick_filename, pdf_filename, png_filename, svg_filename):
    dendrogram_data = {"id": batch_id, "png_file": File(open(png_filename, "rb")),
                       "pdf_file": File(open(pdf_filename, "rb")), "svg_file": File(open(svg_filename, "rb")),
                       "emf_file": File(open(emf_filename, "rb")), "newick_file": File(open(newick_filename, "rb"))}
    serializer = DendrogramSerializer(data=dendrogram_data)
    if serializer.is_valid():
        serializer.save()
    else:
        print(serializer.errors)


@shared_task
def plot_dendrogram(batch_id):
    input_dir = os.path.join(settings.MEDIA_ROOT, "uploads", batch_id)
    output_dir = os.path.join(settings.MEDIA_ROOT, "temp", batch_id)
    files.create_if_not_exist(output_dir)

    emf_filename, newick_filename, pdf_filename, png_filename, svg_filename = plot(input_dir, output_dir)
    save(batch_id, emf_filename, newick_filename, pdf_filename, png_filename, svg_filename)
    shutil.rmtree(output_dir)
