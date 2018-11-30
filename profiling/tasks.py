from __future__ import absolute_import, unicode_literals

import os
import shutil
import zipfile
import pandas as pd
from celery import shared_task
from django.conf import settings
from django.core.files import File

from profiling.serializers import ProfileSerializer
from src.algorithms import profiling
from src.utils import files


def separate_profiles(profile, profile_dir):
    df = pd.read_csv(profile, sep="\t", index_col=0)
    for col in df.columns:
        df[[col]].to_csv(os.path.join(profile_dir, str(col) + ".tsv"), sep="\t")


def zip_folder(filepath):
    filename = filepath + '.zip'
    with zipfile.ZipFile(filename, mode='w', compression=zipfile.ZIP_LZMA) as zip:
        for root, folders, files in os.walk(filepath):
            for f in files:
                zip.write(os.path.join(root, f))
    return filename


def profile(batch_id, database, input_dir, occr_level, output_dir):
    profile_filename = os.path.join(output_dir, "profile.tsv")
    profiling.profiling(output_dir, input_dir, database, occr_level=occr_level, threads=2)

    # Separate individual profiles and zip
    profile_dir = os.path.join(output_dir, batch_id)
    files.create_if_not_exist(profile_dir)
    separate_profiles(profile_filename, profile_dir)
    zip_filename = zip_dir(profile_dir)

    profile_data = {"id": batch_id, "file": File(open(profile_filename, "rb")),
                    "zip": File(open(zip_filename, "rb")),
                    "occurrence": occr_level, "database": database}
    serializer = ProfileSerializer(data=profile_data)
    if serializer.is_valid():
        serializer.save()
    else:
        print(serializer.errors)

    shutil.rmtree(output_dir)
