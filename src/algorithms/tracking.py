import numpy as np
import os
import pandas as pd
import zipfile
from src.utils import nosql


def get_profile(track, biosample):
    return track.find_one({'BioSample': biosample})["profile"]


def distance_against_all(query_profile, track, top_n=100):
    distances = {}
    for subject in track.find({}, {'_id': 0, 'BioSample': 1, 'profile': 1}):
        subject_profile = subject["profile"]
        same_loci = query_profile.keys() & subject_profile.keys()
        diff_alleles_count = sum(map(lambda loci: query_profile[loci] != subject_profile[loci], same_loci))
        diff_loci_count = len((query_profile.keys() | subject_profile.keys()) - same_loci)
        distances[subject['BioSample']] = diff_alleles_count + diff_loci_count
    top_n_dist = sorted(distances.items(), key=lambda x: x[1])[0:top_n]
    return top_n_dist


def save_profiles(samples, track, pf_dir):
    samples = list(map(lambda x: x[0], samples))
    profiles = track.find({'BioSample': {'$in': samples}}, {'_id': 0, 'BioSample': 1, 'profile': 1})
    for profile in profiles:
        s = pd.Series(profile['profile'], name=profile['BioSample'])
        pf = pd.DataFrame(s)
        pf.to_csv(os.path.join(pf_dir, profile['BioSample'] + ".tsv"), sep="\t")


def add_metadata(distances, track):
    samples = list(map(lambda x: x[0], distances))
    metadata = track.find({'BioSample': {'$in': samples}}, {'_id': 0, 'profile': 0})
    metadata = pd.DataFrame(list(metadata))
    distances = pd.DataFrame(distances, columns=['BioSample', 'distances'])
    results = pd.merge(distances, metadata, on='BioSample').sort_values('distances')
    return results


def calculate_distances(profile_filename, track):
    query_profile = pd.read_csv(profile_filename, sep="\t", index_col=0, usecols=[0, 1])
    query_profile = next(query_profile.iteritems())[1].dropna().to_dict()
    distances = distance_against_all(query_profile, track)
    results = add_metadata(distances, track)
    results = results.fillna('')
    return distances, results


def to_json(id, output_dir, results):
    results_file = os.path.join(output_dir, id[0:8] + ".json")
    with open(results_file, "w") as file:
        file.write(results.to_json(orient='records'))
    return results_file


def to_zip(filename, results_file, prof_dir):
    with zipfile.ZipFile(filename, mode='w', compression=zipfile.ZIP_DEFLATED) as zip:
        zip.write(results_file, os.path.basename(results_file))
        for root, folders, files in os.walk(prof_dir):
            for f in files:
                zip.write(os.path.join(root, f),
                          os.path.join(os.path.basename(prof_dir), f))
    return filename


def tracking(id, database, output_dir, profile_filename):
    track = nosql.get_dbtrack(database)
    distances, results = calculate_distances(profile_filename, track)
    results_json = to_json(id, output_dir, results)

    results_file = os.path.join(output_dir, id[0:8] + ".tsv")
    results.to_csv(results_file, sep="\t")

    prof_dir = os.path.join(output_dir, "profiles")
    os.makedirs(prof_dir, exist_ok=True)
    save_profiles(distances, track, prof_dir)

    results_zip = os.path.join(output_dir, id[0:8] + '.zip')
    to_zip(results_zip, results_file, prof_dir)

    return results_json, results_zip
