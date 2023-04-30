#!/usr/bin/env python

import argparse
import urllib.request

BASE_URL = "ftp://ccg.epfl.ch/epdnew/"

parser = argparse.ArgumentParser(description='Fetch EPDNEW data')

parser.add_argument('-g', '--genome', dest='genome', type=str, help='Genome', required=True)
parser.add_argument('-o', '--output', dest='output', type=str, help='Output file', required=True)

args = parser.parse_args()

genomeMap = {
    "hg38": "H_sapiens/006/Hs_EPDnew_006_hg38.bed",
    "hg19": "H_sapiens/006/Hs_EPDnew_006_hg19.bed",
    "mm10": "M_musculus/003/Mm_EPDnew_003_mm10.bed",
    "mm9": "M_musculus/002/Mm_EPDnew_002_mm9.bed",
    "dr7": "D_rerio/001/Dr_EPDnew_001_danRer7.bed",
    "ce6": "C_elegans/001/Ce_EPDnew_001_ce6.bed",
    "dm6": "D_melanogaster/005/Dm_EPDnew_005_dm6.bed",
    "rn6": "R_norvegicus/001/Rn_EPDnew_001_rn6.bed",
    "sacCer3": "S_cerevisiae/002/Sc_EPDnew_002_sacCer3.bed"
}

if args.genome not in genomeMap:
    raise Exception("Genome not supported")

url = BASE_URL + genomeMap[args.genome]

urllib.request.urlretrieve(url, args.output)