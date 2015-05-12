#!/usr/bin/python

import sys
import os
import argparse
import getopt
import pybedtools

#class Error(Exception):
#    """Base-class for exceptions in this module."""

#class UsageError(Error):
#    def __init__(self, msg):
#        self.msg = msg

def main(argv):
    inputfile = ''
    outputfile = ''
    vcf = ''

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--vcf", type=str, dest="vcf", help="Input variant file (vcf)", required=True)
    parser.add_argument("-o", "--output", type=str, dest="out", help="Output file (tabular)", required=True)
    parser.add_argument("-v", "--verbosity", action="count", default=0)

    args = parser.parse_args()

    if args.verbosity >= 2:
        print "{} to the power {} equals {}".format(args.v, args.o, answer)
    elif args.verbosity >= 1:
        print "{}^{} == {}".format(args.x, args.y, answer)
    else:
        print "Starting ..."

    avcf = pybedtools.BedTool(args.vcf)
    current_info = avcf[0][7]
    info = current_info.split("ANN=")
    ann = info[1].split("|")

    for x in xrange(0, len(avcf)):
        current_chr = avcf[x][0]
        current_pos = avcf[x][1]
        current_info = avcf[x][7]
        info_vep = current_info.split("CSQ=")
        info_snpeff = current_info.split("ANN=")

        csq = info_vep[1].split("|")
        current_exon = csq[26]

        ann = info_snpeff[1].split("|")
        current_gene_name = ann[3]
        current_gene2 = ann[4]
        out_str = ["chr"+current_chr, current_pos, ann[1], current_gene_name, current_exon]
        print "\t".join(out_str)

    print len(avcf)


if __name__ == "__main__":
    main(sys.argv[1:])

