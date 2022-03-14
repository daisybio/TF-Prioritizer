import os
import pandas as pd

inputDirectory = "{INPUT_DIRECTORY}"
outputDirectory = "{OUTPUT_DIRECTORY}"
f_outputMeans = "{OUTPUT_MEANS_FILE}"
f_symbolMap = "{SYMBOL_MAP_FILE}"
symbolMap = {}
with open(f_symbolMap, "r") as symbolMapFile:
    for line in symbolMapFile.readlines()[1:]:
        split = line.split("\t")
        symbolMap[split[0]] = split[1]

geneID = "{GENE_ID_FILE}"
df_geneID = pd.read_csv(geneID)
size = df_geneID.size

df_allCombined = df_geneID.copy()
os.makedirs(outputDirectory, exist_ok=True)

for f_sample in os.listdir(inputDirectory):
    absolutePath = os.path.join(inputDirectory, f_sample)
    df_sample = pd.read_csv(absolutePath).astype("int32")

    if not df_sample.size == size:
        exit(1)

    df_combined = df_geneID.join(df_sample)
    f_output = os.path.join(outputDirectory,
                            f_sample.rstrip("txt") + "tsv")
    with open(f_output, "w+") as outputFile:
        df_combined.to_csv(outputFile, sep="\t", index=False, header=True)
    df_allCombined = df_allCombined.join(df_sample)

df_allCombined["MeanCount"] = df_allCombined.mean(axis=1)
df_allCombined["MeanCount"] = df_allCombined["MeanCount"].round(0).astype("int32")
df_allCombined = df_allCombined[["Geneid", "MeanCount"]]
df_allCombined["Symbol"] = df_allCombined["Geneid"].map(symbolMap).fillna("NO_SYMBOL")
columns = df_allCombined.columns.tolist()
df_allCombined = df_allCombined[columns[-1:] + columns[:-1]]
with open(f_outputMeans, "w+") as meansFile:
    df_allCombined.to_csv(meansFile, sep="\t", index=False, header=True)
