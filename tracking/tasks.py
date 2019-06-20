import os
import zipfile
import binascii
import requests
import numpy as np
from celery import shared_task
import pandas as pd
from django.conf import settings
from src.utils import nosql


def get_profile(track, biosample):
    return track.find_one({'BioSample': biosample})["profile"]


def distance_against_all(query_profile, track, top_n=100):
    distances = {}
    for profile in track.find():
        ref_profile = pd.Series(data=profile["profile"])
        df = pd.concat([ref_profile, query_profile], axis=1, keys=['ref_profile', 'query_profile']).fillna('')
        distances[profile['BioSample']] = (df['query_profile'] != df['ref_profile']).sum()
    top_n_dist = pd.Series(distances).sort_values()[0:top_n]
    return top_n_dist


def save_profiles(biosamples, track, prof_dir):
    for biosample in biosamples:
        prof = pd.DataFrame()
        prof[biosample] = pd.Series(data=get_profile(track, biosample))
        prof.to_csv(os.path.join(prof_dir, biosample + ".tsv"), sep="\t")


def add_metadata(distances, track):
    results = pd.DataFrame()
    results['distance'] = distances
    for sample in results.index:
        metadata = track.find_one({'BioSample': sample})
        for col in metadata:
            if col != "profile" and col != "_id":
                results.loc[sample, col] = metadata[col]
    return results


def calculate_distances(profile_filename, track):
    query_profile = pd.read_csv(profile_filename, sep="\t", index_col=0)
    query_profile = query_profile[query_profile.columns[0]]
    distances = distance_against_all(query_profile, track)
    results = add_metadata(distances, track)
    results.replace(np.nan, "", regex=True, inplace=True)
    return distances, results


def calculate_tracking_results(id, database, output_dir, profile_filename):
    track = nosql.get_dbtrack(database)
    distances, results = calculate_distances(profile_filename, track)
    results_json = to_json(id, output_dir, results)

    results_file = os.path.join(output_dir, id[0:8] + ".tsv")
    results.to_csv(results_file, sep="\t")

    prof_dir = os.path.join(output_dir, "profiles")
    os.makedirs(prof_dir, exist_ok=True)
    save_profiles(distances.index, track, prof_dir)

    results_zip = os.path.join(output_dir, id[0:8] + '.zip')
    to_zip(results_zip, results_file, prof_dir)

    return results_json, results_zip


def to_json(id, output_dir, results):
    results_file = os.path.join(output_dir, id[0:8] + ".json")
    with open(results_file, "w") as file:
        file.write(results.to_json(orient='records'))
    return results_file


def to_zip(filename, results_file, prof_dir):
    with zipfile.ZipFile(filename, mode='w', compression=zipfile.ZIP_LZMA) as zip:
        zip.write(results_file, os.path.basename(results_file))
        for root, folders, files in os.walk(prof_dir):
            for f in files:
                zip.write(os.path.join(root, f),
                          os.path.join(os.path.basename(prof_dir), f))
    return filename


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

    json_file, zip_file = calculate_tracking_results(id, database, output_dir, profile_filename)
    save(id, json_file, zip_file, url)
