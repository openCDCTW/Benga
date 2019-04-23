import os
import zipfile
import numpy as np
from collections import Counter
from celery import shared_task
import pandas as pd
from django.conf import settings
from django.core.files import File
from src.utils import nosql
from src.algorithms import profiling
from src.utils import files
from tracking.serializers import TrackedResultsSerializer


def get_profile(track, biosample):
    return track.find_one({'BioSample': biosample})["profile"]


def distance_against_all(query_profile, track, top_n=100):
    distances = pd.Series()
    df = pd.DataFrame()
    df["query_profile"] = query_profile
    for profile in track.find():
        df["ref_profile"] = pd.Series(data=profile["profile"])
        distances.at[profile['BioSample']] = Counter(df["query_profile"] == df["ref_profile"])[0]
    top_n_dist = distances.sort_values()[0:top_n]
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


def to_db(id, results_json, results_zip):
    results = {"id": id, "json": File(open(results_json, "rb")),
               "zip": File(open(results_zip, "rb"))}
    serializer = TrackedResultsSerializer(data=results)
    if serializer.is_valid():
        serializer.save()
    else:
        print(serializer.errors)


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


@shared_task
def track(id, database):
    input_dir = os.path.join(settings.MEDIA_ROOT, "tracking", id)
    profile_filename = os.path.join(input_dir, "profile.tsv")
    output_dir = os.path.join(settings.MEDIA_ROOT, "temp", id)
    os.makedirs(output_dir, exist_ok=True)

    # Prepare tracking result in json format
    query_profile = pd.read_csv(profile_filename, sep="\t", index_col=0)
    query_profile = query_profile[query_profile.columns[0]]
    track = nosql.get_dbtrack(database)
    distances = distance_against_all(query_profile, track)
    results = add_metadata(distances, track)
    results.replace(np.nan, "", regex=True, inplace=True)
    results_json = to_json(id, output_dir, results)

    # Prepare tracking result in tsv format and profiles
    results_file = os.path.join(output_dir, id[0:8] + ".tsv")
    results.to_csv(results_file, sep="\t")
    prof_dir = os.path.join(output_dir, "profiles")
    os.makedirs(prof_dir, exist_ok=True)
    save_profiles(distances.index, track, prof_dir)
    results_zip = os.path.join(output_dir, id[0:8] + '.zip')
    to_zip(results_zip, results_file, prof_dir)

    to_db(id, results_json, results_zip)


@shared_task
def profile_and_track(id, allele_db, occr_level, profile_db):
    input_dir = os.path.join(settings.MEDIA_ROOT, "tracking", id)
    output_dir = os.path.join(settings.MEDIA_ROOT, "temp", id)
    os.makedirs(output_dir, exist_ok=True)

    profile_filename = os.path.join(output_dir, "profile.tsv")
    profiling.profiling(output_dir, input_dir, allele_db, occr_level=occr_level, threads=2)

    query_profile = pd.read_csv(profile_filename, sep="\t", index_col=0)
    query_profile = query_profile[query_profile.columns[0]]
    track = nosql.get_dbtrack(profile_db)
    distances = distance_against_all(query_profile, track)
    results = add_metadata(distances, track)
    results.replace(np.nan, "", regex=True, inplace=True)

    results_file = os.path.join(output_dir, id[0:8] + ".json")
    with open(results_file, "w") as file:
        file.write(results.to_json(orient='records'))
    to_db(id, results_file)
