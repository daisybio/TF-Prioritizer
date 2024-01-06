#!/usr/bin/env python3

import argparse

def convert_bed_to_gff(bed_filename, gff_filename):
    gff_lines = []

    with open(bed_filename, 'r') as bed:
        for line in bed:
            parts = line.strip().split()
            if len(parts) < 3:
                # skip invalid lines
                continue

            seqid, start, end = parts[0], int(parts[1]), parts[2]
            # adjust for 0-based to 1-based coordinate
            start += 1
            
            # Set default values for the other GFF columns (referring to the names in https://en.wikipedia.org/wiki/General_feature_format)
            source = 'bed2gff'
            type = 'region'
            score = '.'
            strand = '.'
            phase = '.'
            attributes = '.'

            gff_line = f"{seqid}\t{source}\t{type}\t{start}\t{end}\t{score}\t{strand}\t{phase}\t{attributes}"
            gff_lines.append(gff_line)

    with open(gff_filename, 'w') as gff:
        gff.write('\n'.join(gff_lines))

def main():
    parser = argparse.ArgumentParser(description='Convert BED file to GFF format.')
    parser.add_argument('--bed', '-i', type=str, help='Input BED file')
    parser.add_argument('--gff', '-o', type=str, help='Output GFF file')

    args = parser.parse_args()

    convert_bed_to_gff(args.bed, args.gff)

if __name__ == "__main__":
    main()
