from time import sleep
import pandas as pd
from pyliftover import LiftOver
import argparse


def convert(c, x, y, s, converter):
    first = converter.convert_coordinate(c, int(x), s)
    second = converter.convert_coordinate(c, int(y), s)

    if (first is None or second is None):
        return None, None

    if len(first) == 0 or len(second) == 0:
        return None, None

    return str(first[0][1]), str(second[0][1])


parser = argparse.ArgumentParser()
parser.add_argument("original", type=argparse.FileType(
    "r"), help="File with positions to uplift from")
parser.add_argument(
    "uplifted", type=argparse.FileType("w"), help="Output file")
parser.add_argument("original_version", type=str)
parser.add_argument("target_version", type=str)

args = parser.parse_args()

df = pd.read_csv(args.original, sep="\t", header=None)

col_chromosome = 3
col_start = 4
col_end = 5
col_strand = 6

converter = LiftOver(args.original_version, args.target_version, search_dir=None, cache_dir=None,
                     use_web=True)

for i, row in df.iterrows():
    chromosome = f"chr{row[col_chromosome]}"
    start_position = row[col_start]
    end_position = row[col_end]
    strand = row[col_strand]

    chromosome = "chrM" if chromosome == "chrMT" else chromosome

    start_converted, end_converted = convert(
        chromosome, start_position, end_position, strand, converter)

    if (start_converted is None or end_converted is None):
        continue

    df.loc[i, col_start] = start_converted
    df.loc[i, col_end] = end_converted

df.to_csv(args.uplifted, sep="\t", header=False, index=False)
