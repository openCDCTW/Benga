def generate_encoder(alleles):
    encoder = {}
    for i, x in enumerate(sorted(set(alleles)), 1):
        if x == '0':
            encoder[x] = 0
        else:
            encoder[x] = i
    return encoder


def to_bionumerics_format(profile):
    wgmlst = profile.fillna('0').transpose()
    for col in wgmlst.columns:
        alleles = wgmlst[col]
        encoder = generate_encoder(alleles)
        wgmlst[col] = [encoder[i] for i in alleles]
    wgmlst.insert(0, 'Key', wgmlst.index)
    return wgmlst

