#!/usr/bin/python

import sys
import getopt
import os
import pybedtools


def main(argv):
    inputfile = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
    except getopt.GetoptError:
        print 'readvcf.py -i <inputfile> -o <outputfile>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'readvcf.py -i <inputfile> -o <outputfile>'
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
    print 'Input VCF file is "', inputfile
    print 'Output file is "', outputfile

if __name__ == "__main__":
   main(sys.argv[1:])


cpath = os.path.dirname(os.path.abspath(__file__))

#kgbed = pybedtools.BedTool('/Users/khalidm/Development/repositories/VarPub/data/1kg.phase1.snp.bed.gz')
kgbed = pybedtools.BedTool('/Users/khalidm/Development/repositories/VarPub/data/1kg.phase1.snp.bed.gz')
#avcf = pybedtools.BedTool('/Users/khalidm/Development/repositories/VarPub/data/BRCA2.snpeff.vcf')
avcf = pybedtools.BedTool('/Users/khalidm/Development/repositories/VarPub/data/BRCA2.snpeff.vcf')

# print "cat a.bed\n" + str(avcf)
# print "cat b.bed\n" + str(bfile)

# print intersect
# print avcf.intersect(kgbed)
# print length
# print len(a)

print avcf[0][7] + " -> " + avcf[0][8]

