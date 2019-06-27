from __future__ import absolute_import, unicode_literals

import os
import shutil
import zipfile
import binascii
import requests
import pandas as pd
from celery import shared_task
from django.conf import settings
from src.algorithms import profiling


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


def save(batch_id, database, occr_level, zip_filename, url):
    profile_data = {"id": batch_id, "occurrence": occr_level, "database": database}
    files = {"zip": open(zip_filename, "rb")}
    r = requests.post(url, data=profile_data, files=files)
    if r.status_code != 201:
        print(r.status_code)


def get_file_number(dir, ext=".tsv"):
    return len(list(filter(lambda x: x.endswith(ext), os.listdir(dir))))


@shared_task
def single_profiling(seq_id, batch_id, database, occr_level, seq_num, filename, sequence, url):
    input_dir = os.path.join(settings.CELERY_ROOT, "uploads", batch_id, seq_id)
    os.makedirs(input_dir, exist_ok=True)
    with open(os.path.join(input_dir, filename), "wb") as file:
        file.write(binascii.a2b_base64(sequence))

    output_dir = os.path.join(settings.CELERY_ROOT, "temp", batch_id)
    os.makedirs(output_dir, exist_ok=True)
    profiling.profiling(output_dir, input_dir, database, occr_level=occr_level,
                        threads=2, profile_file=filename, generate_bn=False)
    return batch_id, output_dir, database, occr_level, seq_num, url


@shared_task
def zip_save(args):
    batch_id, output_dir, database, occr_level, seq_num, url = args
    if os.path.exists(output_dir) and seq_num == get_file_number(output_dir):
        zip_filename = zip_folder(output_dir)
        save(batch_id, database, occr_level, zip_filename, url)
        shutil.rmtree(output_dir)

