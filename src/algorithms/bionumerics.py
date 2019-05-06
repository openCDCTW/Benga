import pandas as pd


def to_bionumerics_format(profile):
    profile = profile.fillna(0)
    profile_values = profile.values
    new_values = []
    for value in profile_values:
        allele_set = set(value)
        allele_set.discard(0)
        encoder = {allele: index for index, allele in enumerate(allele_set, 1)}
        value = [encoder.get(i, 0) for i in value]
        new_values.append(value)
    df = pd.DataFrame(data=new_values, index=profile.index, columns=profile.columns)
    df = df.transpose()
    df.insert(0, 'Key', df.index)
    return df
