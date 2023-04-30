#!/usr/bin/env python

import json
import argparse

parser = argparse.ArgumentParser(description='Convert TF groups to single TFs')

parser.add_argument('-i', '--input', dest='input', type=argparse.FileType('r'), help='Input file', required=True)
parser.add_argument('-o', '--output', dest='output', type=argparse.FileType('w'), help='Output file', required=True)
parser.add_argument('-j', '--json', dest='json', type=argparse.FileType('w'), help='JSON file', required=True)

args = parser.parse_args()

tf_groups = args.input.read().splitlines()

dictionary = { group: group.split('..') for group in tf_groups }
tfs = list(set([ tf for tfs in dictionary.values() for tf in tfs ]))
tfs.sort()

json.dump(dictionary, args.json, indent=4)

for tf in tfs:
    args.output.write(tf + '\n')