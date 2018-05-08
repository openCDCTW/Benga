import pandas as pd

def bionumerics(profile):
    wgmlst= profile.fillna('0').transpose()
    for col in wgmlst.columns:
        allele_id_pairs = {}
        c = 1
        for i in sorted(set(wgmlst[col])):
            if i == '0':
                allele_id_pairs[i] = 0
            else:
                allele_id_pairs[i] = c
                c += 1
        wgmlst[col] = [allele_id_pairs[i] for i in wgmlst[col]]
    wgmlst.insert(0, 'Key', wgmlst.index)
    return wgmlst