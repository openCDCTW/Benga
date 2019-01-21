import os
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


def distance_against_all(query_profile, track, top_n=100):
    distances = pd.Series()
    df = pd.DataFrame()
    df["query_profile"] = query_profile
    for profile in track.find():
        df["ref_profile"] = pd.Series(data=profile["profile"])
        distances.at[profile['BioSample']] = Counter(df["query_profile"] == df["ref_profile"])[0]
    return distances.sort_values()[0:top_n]


def add_metadata(distances, track):
    results = pd.DataFrame()
    results['distance'] = distances
    for sample in results.index:
        metadata = track.find_one({'BioSample': sample})
        for col in metadata:
            if col != "profile" and col != "_id":
                results.loc[sample, col] = metadata[col]
    return results


def to_db(id, results_file):
    results = {"id": id, "json": File(open(results_file, "rb"))}
    serializer = TrackedResultsSerializer(data=results)
    if serializer.is_valid():
        serializer.save()
    else:
        print(serializer.errors)


@shared_task
def track(query_profile, database):
    track = nosql.get_dbtrack(database)
    distances = distance_against_all(query_profile, track)
    results = add_metadata(distances, track)
    return results


@shared_task
def profile_and_track(id, allele_db, occr_level, profile_db):
    input_dir = os.path.join(settings.MEDIA_ROOT, "tracking", id)
    output_dir = os.path.join(settings.MEDIA_ROOT, "temp", id)
    files.create_if_not_exist(output_dir)

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
