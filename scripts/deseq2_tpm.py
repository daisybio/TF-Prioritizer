import pandas as pd

countFileName = "{COUNTFILE}"
lengthFileName = "{LENGTHSFILE}"
targetFileName = "{TARGETFILE}"

df = pd.read_csv(countFileName, sep="\t")
mapping = {}

with open(lengthFileName, "r") as lengthFile:
    lines = lengthFile.readlines()
    ensg_index = lines[0].split("\t").index("ENSG")
    length_index = lines[0].split("\t").index("length")

    for entry in lines[1:]:
        split = entry.split("\t")
        mapping[split[ensg_index]] = split[length_index]

df["length"] = df["Geneid"].map(mapping)
df = df[df["length"] != "NA"]
df["length"] = df["length"].astype(int)

df["rpk"] = df["MeanCount"] / df["length"]

rpk_sum_per_million = sum(df["rpk"]) / 1e6

df["tpm"] = df["rpk"] / rpk_sum_per_million

remaining_columns = df.columns.tolist()
remaining_columns.remove("rpk")
df = df[remaining_columns]

df = df.sort_values("Geneid")

with open(targetFileName, "w+") as targetFile:
    df.to_csv(targetFile, sep="\t", index=False)
