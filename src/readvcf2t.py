'''
A tool to annotate and print variants in tabular format.
Author: Khalid Mahmood (khalid.mahmood@unimelb.edu.au).
Copyright: 2015
'''

#!/usr/bin/python


from utils import findlist

import sys
import os
import argparse
import getopt
import vcf
import array
import pysam

#class Error(Exception):
#    """Base-class for exceptions in this module."""

#class UsageError(Error):
#    def __init__(self, msg):
#        self.msg = msg

# MAIN

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
    cadd_tbx = pysam.TabixFile("data/whole_genome_SNVs_inclAnno.tsv.gz")
    cadd_indel_tbx = pysam.TabixFile("data/InDels_inclAnno.tsv.gz")
    fathmm_tbx = pysam.TabixFile("data/fathmm-MKL_Current_zerobased.tab.gz")
    exac_tbx = pysam.TabixFile("data/ExAC.r0.3.sites.vep.vcf.gz")

    outputfile.write("chr\tpos\tid\tref\talt\tannotation\tgene_name\tlof" \
            "\texon\taa_pos\tpoly/sift\tAF\tGMAF\t1kgEMAF\tESPEMAF\t" \
            #"HETEUR\tHOMEUR\t
            "ExAC_AF\tExAC_EAS\tExAC_NFE\tExAC_FIN\tExAC_SAS\tExAC_AFR\tExAC_AMR\tExAC_OTH\t" \
            "CADD\tmaxCADD\tpriPhCons\tGerpRS\t" \
            "FATHMM\n")

    vcf_reader = vcf.Reader(open(args.vcf, 'r'))

    for record in vcf_reader:
        outputfile.write("\t".join(annotator(record)))
        outputfile.write("\n")

    outputfile.close()


if __name__ == "__main__":
    main(sys.argv)


