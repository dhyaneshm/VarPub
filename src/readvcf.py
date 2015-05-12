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

    #answer = args.x**args.y
    #print "INPUT = " , args.vcf

    if args.verbosity >= 2:
        print "{} to the power {} equals {}".format(args.v, args.o, answer)
    elif args.verbosity >= 1:
        print "{}^{} == {}".format(args.x, args.y, answer)
    else:
        print "Starting ..."

    avcf = pybedtools.BedTool(args.vcf)
    print avcf[0][7] + " -> " + avcf[0][8]

if __name__ == "__main__":
    main(sys.argv[1:])

