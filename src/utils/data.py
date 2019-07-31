import pandas as pd


def convert_data_type(profile):
    profile = profile.fillna(0)
    data = []
    for str_value in profile.values:
        str_value_set = set(str_value)
        str_value_set.discard(0)
        encoder = {string: _ for _, string in enumerate(str_value_set, 1)}
        int_value = [encoder.get(i, 0) for i in str_value]
        data.append(int_value)
    df = pd.DataFrame(data=data, index=profile.index, columns=profile.columns)
    return df
