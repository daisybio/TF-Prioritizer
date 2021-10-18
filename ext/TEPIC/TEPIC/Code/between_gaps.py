#not needed anymore
#this is a corpse
import sys

def read_regions(infile):
    peak_lines=[]
    for l in infile:
        peak_lines.append(l)

    return peak_lines

def main():
    infile=open(sys.argv[1],"r")
    peaks=read_regions(infile)

    outfile=sys.argv[2]

main()
