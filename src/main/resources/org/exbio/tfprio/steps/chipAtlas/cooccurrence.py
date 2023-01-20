import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('input', type=str)
parser.add_argument('output', type=str)

args = parser.parse_args()

df = pd.read_table(args.input, header=None, usecols=[3], names=["TF"])
df = df["TF"].str.get_dummies(sep="|")

matrix = df.T.dot(df)

matrix.to_csv(args.output, sep="\t")
