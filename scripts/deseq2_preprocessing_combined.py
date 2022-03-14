import os
import pandas as pd

inputDirectory1 = "{INPUT_DIRECTORY1}"
inputDirectory2 = "{INPUT_DIRECTORY2}"
outputFile = "{OUTPUT_FILE}"
geneIdFile = "{GENE_ID}"

df_combined = pd.read_csv(geneIdFile)

for fileName_sample in os.listdir(inputDirectory1):
    absolutePath = os.path.join(inputDirectory1, fileName_sample)
    df_sample = pd.read_csv(absolutePath, sep="\t")

    if df_sample["Geneid"].equals(df_combined["Geneid"]):
        columns = df_sample.columns.tolist()
        columns.remove("Geneid")
        df_combined = df_combined.join(df_sample[columns])

for fileName_sample in os.listdir(inputDirectory2):
    absolutePath = os.path.join(inputDirectory2, fileName_sample)
    df_sample = pd.read_csv(absolutePath, sep="\t")

    if df_sample["Geneid"].equals(df_combined["Geneid"]):
        columns = df_sample.columns.tolist()
        columns.remove("Geneid")
        df_combined = df_combined.join(df_sample[columns])

with open(outputFile, "w+") as file:
    df_combined.to_csv(file, index=False, sep="\t")
