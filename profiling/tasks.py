from __future__ import absolute_import, unicode_literals

import os
import shutil
import zipfile
import pandas as pd
from celery import shared_task
from django.conf import settings
from django.core.files import File

from profiling.serializers import ProfileSerializer
from profiling.models import Batch, Sequence
from src.algorithms import profiling
import dendrogram.tasks as tree


def separate_profiles(profile, profile_dir):
    df = pd.read_csv(profile, sep="\t", index_col=0)
    for col in df.columns:
        df[[col]].to_csv(os.path.join(profile_dir, str(col) + ".tsv"), sep="\t")


def zip_folder(filepath):
    filename = filepath + '.zip'
    with zipfile.ZipFile(filename, mode='w', compression=zipfile.ZIP_LZMA) as zip:
        for root, folders, files in os.walk(filepath):
            for f in files:
                zip.write(os.path.join(root, f), f)
    return filename


def profile(batch_id, database, input_dir, occr_level, output_dir):
    profile_filename = os.path.join(output_dir, "profile.tsv")
    profiling.profiling(output_dir, input_dir, database, occr_level=occr_level, threads=2)

    # Separate individual profiles and zip
    profile_dir = os.path.join(output_dir, batch_id)
    os.makedirs(profile_dir, exist_ok=True)
    separate_profiles(profile_filename, profile_dir)
    zip_filename = zip_folder(profile_dir)
    return profile_filename, zip_filename


def save(batch_id, database, occr_level, zip_filename):
    profile_data = {"id": batch_id, "zip": File(open(zip_filename, "rb")),
                    "occurrence": occr_level, "database": database}
    serializer = ProfileSerializer(data=profile_data)
    if serializer.is_valid():
        serializer.save()


def get_filename(seq_id):
    return os.path.basename(Sequence.objects.get(pk=seq_id).file.name)


def get_seq_number(batch_id):
    return Batch.objects.get(pk=batch_id).seq_num


def get_file_number(dir, ext=".tsv"):
    return len(list(filter(lambda x: x.endswith(ext), os.listdir(dir))))


@shared_task
def single_profiling(seq_id, batch_id, database, occr_level):
    input_dir = os.path.join(settings.MEDIA_ROOT, "uploads", batch_id, seq_id)
    output_dir = os.path.join(settings.MEDIA_ROOT, "temp", batch_id)
    os.makedirs(output_dir, exist_ok=True)
    profiling.profiling(output_dir, input_dir, database, occr_level=occr_level,
                        threads=2, profile_file=get_filename(seq_id))
    return batch_id, output_dir, database, occr_level


@shared_task
def zip_save(args):
    batch_id, output_dir, database, occr_level = args
    seq_num = get_seq_number(batch_id)
    if seq_num == get_file_number(output_dir):
        zip_filename = zip_folder(output_dir)
        save(batch_id, database, occr_level, zip_filename)
        shutil.rmtree(output_dir)


@shared_task
def batch_profiling(batch_id, database, occr_level):
    input_dir = os.path.join(settings.MEDIA_ROOT, "uploads", batch_id)
    output_dir = os.path.join(settings.MEDIA_ROOT, "temp", batch_id)
    os.makedirs(output_dir, exist_ok=True)

    _, zip_filename = profile(batch_id, database, input_dir, occr_level, output_dir)
    save(batch_id, database, occr_level, zip_filename)
    shutil.rmtree(output_dir)


@shared_task
def profile_and_tree(batch_id, database, occr_level):
    input_dir = os.path.join(settings.MEDIA_ROOT, "uploads", batch_id)
    output_dir = os.path.join(settings.MEDIA_ROOT, "temp", batch_id)
    os.makedirs(output_dir, exist_ok=True)

    # profile
    profile_filename, zip_filename = profile(batch_id, database, input_dir, occr_level, output_dir)
    save(batch_id, database, occr_level, profile_filename, zip_filename)

    # plot dendrogram
    emf_filename, newick_filename, pdf_filename, png_filename, svg_filename = tree.plot(output_dir, output_dir)
    tree.save(batch_id, emf_filename, newick_filename, pdf_filename, png_filename, svg_filename)

    shutil.rmtree(output_dir)
