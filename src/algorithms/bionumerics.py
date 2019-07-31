from ..utils.data import convert_data_type


def to_bionumerics_format(profile):
    profile = integer_encoding(profile).T
    profile.insert(0, 'Key', profile.index)
    return profile
