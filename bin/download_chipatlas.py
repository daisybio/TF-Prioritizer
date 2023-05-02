#!/usr/bin/env python

import argparse
import pandas as pd
import requests
from io import StringIO

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--output', required=True, help='Output file name')
parser.add_argument('-u', '--url', required=True, help='URL to download')

args = parser.parse_args()

text = requests.get(args.url).text.splitlines()[1:]
text = [line.split('\t') for line in text]

if not text:
    print('No data found')
    exit(0)

df = pd.DataFrame(text)
df[3] = df[3].apply(lambda x: x[0:x.find(';')])

df.to_csv(args.output, sep='\t', header=None, index=False)