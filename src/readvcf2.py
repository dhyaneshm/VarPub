#!/usr/bin/python

import sys
import os
import argparse
import getopt
#import pybedtools
import vcf
import array

#class Error(Exception):
#    """Base-class for exceptions in this module."""

#class UsageError(Error):
#    def __init__(self, msg):
#        self.msg = msg

def main(argv):

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--vcf", type=str, dest="vcf", help="Input variant file (vcf)", required=True)
    parser.add_argument("-o", "--output", type=str, dest="out", help="Output file (tabular)", required=True)
    parser.add_argument("-v", "--verbosity", action="count", default=0)

    args = parser.parse_args()
    outputfile = open(args.out, "w")

    if args.verbosity >= 2:
        print "{} to the power {} equals {}".format(args.v, args.o, answer)
    elif args.verbosity >= 1:
        print "{}^{} == {}".format(args.x, args.y, answer)
    #else:
    #    print "Starting ..."

    outputfile.write("chr\tpos\tref\talt\tannotation\tgene_name\tlof\texon\taa_pos\tpoly/sift\teur_maf")


    vcf_reader = vcf.Reader(open(args.vcf, 'r'))
    for record in vcf_reader:
        current_chr = record.CHROM
        current_pos = record.POS
        current_ref = record.REF
        current_alt = ','.join(str(v) for v in record.ALT)

        # VEP
        if "CSQ" in record.INFO:
            current_vep = record.INFO['CSQ'][0].split('|')
            current_feature = current_vep[2]
            current_feature_type = current_vep[3]
            current_consequence = current_vep[4]
        else:
            current_feature = ''
            current_feature_type = ''
            current_consequence = ''


        # SnpEff
        current_ann = record.INFO['ANN'][0].split('|')
        current_annotation = current_ann[1]
        # current_annotation_impact = current_ann[2]
        current_gene = current_ann[3]


        out_str = ["chr"+current_chr, str(current_pos), current_ref, current_alt,
                current_annotation, current_gene ]
        out_str = [x or '.' for x in out_str]
        outputfile.write("\t".join(out_str))
        outputfile.write("\n")

    outputfile.close()



if __name__ == "__main__":
    main(sys.argv)

