from __future__ import absolute_import, unicode_literals
from celery import shared_task
import os
from django.conf import settings
import shutil
from django.core.files import File
from src.algorithms import profiling
from src.utils import files
from profiling.serializers import ProfileSerializer


@shared_task
def do_profiling(batch_id, database, occr_level):
    input_dir = os.path.join(settings.MEDIA_ROOT, "uploaded_sequences", batch_id)
    output_dir = os.path.join(settings.MEDIA_ROOT, "temp", batch_id)
    files.create_if_not_exist(output_dir)

    profile_filename = os.path.join(output_dir, "profile.tsv")
    profiling.profiling(output_dir, input_dir, database, occr_level=occr_level, threads=2)

    profile_data = {"id": batch_id, "file": File(open(profile_filename, "rb")),
                    "occurrence": occr_level, "database": database}
    serializer = ProfileSerializer(data=profile_data)
    if serializer.is_valid():
        serializer.save()
    else:
        print(serializer.errors)

    shutil.rmtree(output_dir)
