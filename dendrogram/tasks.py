from __future__ import absolute_import, unicode_literals

import os
import pandas as pd
import shutil
import subprocess
from celery import shared_task
from django.conf import settings
from django.core.files import File

from dendrogram.serializers import DendrogramSerializer
from dendrogram.models import Batch
from src.algorithms import clustering


def read_profiles(input_dir):
    files = list(filter(lambda x: x.endswith(".tsv"), os.listdir(input_dir)))
    if len(files) == 1:
        return pd.read_csv(os.path.join(input_dir, files[0]), sep="\t", index_col=0)
    else:
        profiles = []
        for filename in files:
            p = pd.read_csv(os.path.join(input_dir, filename), sep="\t", index_col=0)
            profiles.append(p)
        return pd.concat(profiles, axis=1, join='inner', sort=False)


def plot(input_dir, output_dir, linkage):
    profiles = read_profiles(input_dir)
    dm = clustering.DistanceMatrix(profiles)
    dendro = clustering.Dendrogram(dm, linkage)
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


def save(batch_id, linkage, emf_filename, newick_filename, pdf_filename, png_filename, svg_filename):
    dendrogram_data = {"id": batch_id, "linkage": linkage, "png_file": File(open(png_filename, "rb")),
                       "pdf_file": File(open(pdf_filename, "rb")), "svg_file": File(open(svg_filename, "rb")),
                       "emf_file": File(open(emf_filename, "rb")), "newick_file": File(open(newick_filename, "rb"))}
    serializer = DendrogramSerializer(data=dendrogram_data)
    if serializer.is_valid():
        serializer.save()
    else:
        print(serializer.errors)


def get_prof_number(batch_id):
    return Batch.objects.get(pk=batch_id).prof_num


def get_linkage(batch_id):
    return Batch.objects.get(pk=batch_id).linkage


def get_file_number(dir, ext=".tsv"):
    return len(list(filter(lambda x: x.endswith(ext), os.listdir(dir))))


@shared_task
def plot_dendrogram(batch_id):
    input_dir = os.path.join(settings.MEDIA_ROOT, "uploads", batch_id)
    output_dir = os.path.join(settings.MEDIA_ROOT, "temp", batch_id)
    os.makedirs(output_dir, exist_ok=True)
    prof_num = get_prof_number(batch_id)
    if prof_num == get_file_number(input_dir):
        linkage = get_linkage(batch_id)
        emf_file, newick_file, pdf_file, png_file, svg_file = plot(input_dir, output_dir, linkage)
        save(batch_id, linkage, emf_file, newick_file, pdf_file, png_file, svg_file)
        shutil.rmtree(output_dir)
