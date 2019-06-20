from __future__ import absolute_import, unicode_literals

import os
import pandas as pd
import shutil
import binascii
import requests
from celery import shared_task
from django.conf import settings
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
        return pd.concat(profiles, axis=1, sort=False)


def plot(input_dir, output_dir, linkage):
    profiles = read_profiles(input_dir)
    dm = clustering.DistanceMatrix(profiles)
    dendro = clustering.Dendrogram(dm, linkage)
    filenames = {}
    filenames["newick"] = os.path.join(output_dir, "dendrogram.newick")
    dendro.to_newick(filenames["newick"])
    filenames["pdf"] = os.path.join(output_dir, "dendrogram.pdf")
    dendro.scipy_tree(filenames["pdf"])
    filenames["svg"] = os.path.join(output_dir, "dendrogram.svg")
    dendro.scipy_tree(filenames["svg"])
    filenames["png"] = os.path.join(output_dir, "dendrogram.png")
    dendro.scipy_tree(filenames["png"])
    return filenames


def save(batch_id, linkage, filenames, url):
    dendrogram_data = {"id": batch_id, "linkage": linkage}
    files = {"png_file": open(filenames["png"], "rb"), "pdf_file": open(filenames["pdf"], "rb"),
             "svg_file": open(filenames["svg"], "rb"), "newick_file": open(filenames["newick"], "rb")}
    r = requests.post(url, data=dendrogram_data, files=files)
    if r.status_code != 201:
        print(r.status_code)


def get_file_number(dir, ext=".tsv"):
    return len(list(filter(lambda x: x.endswith(ext), os.listdir(dir))))


@shared_task
def plot_dendrogram(batch_id, linkage, prof_num, profile, filename, url):
    input_dir = os.path.join(settings.CELERY_ROOT, "uploads", batch_id)
    os.makedirs(input_dir, exist_ok=True)
    with open(os.path.join(input_dir, filename), "wb") as file:
        file.write(binascii.a2b_base64(profile))

    output_dir = os.path.join(settings.CELERY_ROOT, "temp", batch_id)
    os.makedirs(output_dir, exist_ok=True)
    if prof_num == get_file_number(input_dir):
        filenames = plot(input_dir, output_dir, linkage)
        save(batch_id, linkage, filenames, url)
        shutil.rmtree(output_dir)
