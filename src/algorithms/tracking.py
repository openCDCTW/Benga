import os
import pandas as pd
import zipfile
from src.utils import nosql


def genetic_distance(query, subject):
    loci = set(query) | set(subject)
    distance = sum(query.get(locus_id) != subject.get(locus_id) for locus_id in loci)
    return distance


def save_profiles(profiles, pf_dir):
    for profile in profiles:
        outfile = os.path.join(pf_dir, profile.name + ".tsv")
        profile.to_csv(outfile, sep="\t", header=True)


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
    nosql_client = nosql.NoSQL(database_name=database)
    core_genome = nosql_client.core_genome

    profile = pd.read_csv(profile_filename, sep='\t', index_col=0, usecols=[0, 1])
    query = profile.iloc[:, 0].dropna().to_dict()
    distances = {doc["NCBIAccession"]: genetic_distance(query, doc["profile"])
                 for doc in nosql_client.profilesIterator()}
    top_n_dist = sorted(distances.items(), key=lambda x: x[1])[:100]
    accs = list(map(lambda x: x[0], top_n_dist))
    metadata = pd.DataFrame(nosql_client.fetch_attrib({'NCBIAccession': {'$in': accs}}, {'_id': 0, 'profile': 0}))
    metadata["Distance(Loci)"] = metadata['NCBIAccession'].map(distances)
    metadata = metadata.sort_values("Distance(Loci)")

    mismatch = dict()
    profiles = []
    for doc in nosql_client.fetch_attrib({'NCBIAccession': {'$in': accs}}, {'_id': 0, 'profile': 1, 'NCBIAccession': 1}):
        mismatch[doc['NCBIAccession']] = len(set(core_genome) - set(doc['profile']))
        profiles.append(pd.Series(doc['profile'], name=doc['NCBIAccession']))
    metadata["VoidLoci"] = metadata['NCBIAccession'].map(mismatch)
    results = metadata.reindex(nosql_client.fields, axis=1).fillna("")

    nosql_client.disconnect()

    results_json = to_json(id, output_dir, results)

    results_file = os.path.join(output_dir, id[0:8] + ".tsv")
    results.to_csv(results_file, sep="\t", index=False)

    pf_dir = os.path.join(output_dir, "profiles")
    os.makedirs(pf_dir, exist_ok=True)
    save_profiles(profiles, pf_dir)

    results_zip = os.path.join(output_dir, id[0:8] + '.zip')
    to_zip(results_zip, results_file, pf_dir)

    return results_json, results_zip
