def integer_encoding(values):
    encoder = {j: i for i, j in enumerate(values.dropna().unique(), 1)}
    return values.map(encoder).fillna(0)
