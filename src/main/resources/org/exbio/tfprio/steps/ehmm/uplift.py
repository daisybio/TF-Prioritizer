import argparse
import pandas as pd
from pyliftover import LiftOver


def convert(converter, chromosome, start, end, strand=None):
    if strand is None:
        first = converter.convert_coordinate(chromosome, int(start))
        second = converter.convert_coordinate(chromosome, int(end))
    else:
        first = converter.convert_coordinate(chromosome, int(start), strand)
        second = converter.convert_coordinate(chromosome, int(end), strand)

    if first is None or second is None:
        return None, None

    if len(first) == 0 or len(second) == 0:
        return None, None

    if first[0][1] - second[0][1] > 0:
        return None, None

    return str(first[0][1]), str(second[0][1])


parser = argparse.ArgumentParser()
parser.add_argument("original", type=argparse.FileType(
    "r"), help="File with positions to uplift from")
parser.add_argument(
    "uplifted", type=argparse.FileType("w"), help="Output file")
parser.add_argument("original_version", type=str)
parser.add_argument("target_version", type=str)
parser.add_argument("col_chromosome", type=int)
parser.add_argument("col_start", type=int)
parser.add_argument("col_end", type=int)

args = parser.parse_args()

df = pd.read_csv(args.original, sep="\t", header=None)

col_chromosome = args.col_chromosome
col_start = args.col_start
col_end = args.col_end

converter = LiftOver(args.original_version, args.target_version, search_dir=None, cache_dir=None,
                     use_web=True)

for i, row in df.iterrows():
    chromosome = str(row[col_chromosome])
    chromosome = f"chr{chromosome}" if not chromosome.startswith("chr") else chromosome
    start_position = row[col_start]
    end_position = row[col_end]

    chromosome = "chrM" if chromosome == "chrMT" else chromosome

    start_converted, end_converted = convert(converter,
                                             chromosome, start_position, end_position)

    if start_converted is None or end_converted is None:
        continue

    df.loc[i, col_start] = start_converted
    df.loc[i, col_end] = end_converted

df.to_csv(args.uplifted, sep="\t", header=False, index=False)
