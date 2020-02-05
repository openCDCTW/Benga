from ..utils.data import integer_encoding


def to_bionumerics_format(profile):
    profile = profile.apply(integer_encoding, axis=1)
    profile.insert(0, 'Key', profile.index)
    return profile
