import os
from celery import shared_task
from django.conf import settings
from django.core.files import File
from src.algorithms.tracking import tracking
from tracking.serializers import TrackedResultsSerializer


def save(id, results_json, results_zip):
    results = {"id": id, "json": File(open(results_json, "rb")),
               "zip": File(open(results_zip, "rb"))}
    serializer = TrackedResultsSerializer(data=results)
    if serializer.is_valid():
        serializer.save()
    else:
        print(serializer.errors)


@shared_task
def track(id, database):
    input_dir = os.path.join(settings.MEDIA_ROOT, "tracking", id)
    profile_filename = os.path.join(input_dir, "profile.tsv")
    output_dir = os.path.join(settings.MEDIA_ROOT, "temp", id)
    os.makedirs(output_dir, exist_ok=True)

    json_file, zip_file = tracking(id, database, output_dir, profile_filename)
    save(id, json_file, zip_file)

