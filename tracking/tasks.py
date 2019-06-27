import os
import binascii
import requests
from celery import shared_task
from django.conf import settings
from src.algorithms.tracking import tracking


def save(id, results_json, results_zip, url):
    tracking_data = {"id": id}
    files = {"json": open(results_json, "rb"), "zip": open(results_zip, "rb")}
    r = requests.post(url, data=tracking_data, files=files)
    if r.status_code != 201:
        print(r.status_code)


@shared_task
def track(id, database, profile, url):
    input_dir = os.path.join(settings.CELERY_ROOT, "tracking", id)
    profile_filename = os.path.join(input_dir, "profile.tsv")
    with open(profile_filename, "wb") as file:
        file.write(binascii.a2b_base64(profile))

    output_dir = os.path.join(settings.CELERY_ROOT, "temp", id)
    os.makedirs(output_dir, exist_ok=True)

    json_file, zip_file = tracking(id, database, output_dir, profile_filename)
    save(id, json_file, zip_file, url)
