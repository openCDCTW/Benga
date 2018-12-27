import os
import json
from collections import Counter
from celery import shared_task
import pandas as pd
from django.conf import settings
from django.core.files import File
from src.utils import nosql
from src.algorithms import profiling
from src.utils import files
from tracking.serializers import TrackedResultsSerializer


def distance_against_all(query_profile, track):
    distances = pd.Series()
    for profile in track.find():
        ref_profile = pd.Series(data=profile["profile"])
        distances.at[profile['BioSample']] = Counter(query_profile == ref_profile)[0]
    return distances


def add_metadata(distances, track, top_n=100):
    results = pd.DataFrame()
    results['distance'] = distances.sort_values()[0:top_n]
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

    results_file = os.path.join(output_dir, id + ".json")
    with open(results_file, "w") as file:
        json_content = json.dumps(results.to_dict('records'))
        file.write(json_content)
    to_db(id, results_file)
