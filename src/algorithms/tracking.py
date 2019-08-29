import numpy as np
import os
import pandas as pd
import zipfile
from src.utils import nosql


def get_profile(track, biosample):
    return track.find_one({'BioSample': biosample})["profile"]


def distance_against_all(query_profile, track, top_n=100):
    query_alleles = query_profile.dropna().to_dict()
    distances = {}
    for profile in track.find():
        ref_alleles = profile["profile"]
        same_loci = query_alleles.keys() & ref_alleles.keys()
        diff_alleles_count = sum(True if query_alleles[loci] != ref_alleles[loci] else False for loci in same_loci)
        diff_loci_count = len((query_alleles.keys()|ref_alleles.keys()) - same_loci)
        distances[profile['BioSample']] = diff_alleles_count + diff_loci_count
    top_n_dist = pd.Series(distances, name='distance').sort_values()[0:top_n]
    return top_n_dist


def save_profiles(biosamples, track, prof_dir):
    for biosample in biosamples:
        prof = pd.DataFrame()
        prof[biosample] = pd.Series(data=get_profile(track, biosample))
        prof.to_csv(os.path.join(prof_dir, biosample + ".tsv"), sep="\t")


def add_metadata(distances, track):
    metadata_set = []
    for sample in distances.index:
        metadata = track.find_one({'BioSample': sample}, {'_id': 0, 'profile': 0})
        metadata_set.append(metadata)
    metadata_table = pd.DataFrame(metadata_set)
    results = pd.merge(distances, metadata_table, left_index=True, right_on='BioSample')
    return results


def calculate_distances(profile_filename, track):
    query_profile = pd.read_csv(profile_filename, sep="\t", index_col=0)
    query_profile = query_profile[query_profile.columns[0]]
    distances = distance_against_all(query_profile, track)
    results = add_metadata(distances, track)
    results.replace(np.nan, "", regex=True, inplace=True)
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
    save_profiles(distances.index, track, prof_dir)

    results_zip = os.path.join(output_dir, id[0:8] + '.zip')
    to_zip(results_zip, results_file, prof_dir)

    return results_json, results_zip
