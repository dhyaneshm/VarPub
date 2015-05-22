#!/usr/bin/python

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

def getcadd(cadd_tbx, current_chr, current_pos):
    for row in cadd_tbx.fetch(current_chr, current_pos-1, current_pos):
        row_info = row.split("\t")
        cadd_ref = row_info[2]
        cadd_alt = row_info[4]
        if(cadd_ref == current_ref and cadd_alt == current_alt):
            cadd_phred = row_info[115]
            cadd_priPhCons = row_info[18]
            cadd_GerpRS = row_info[26]
        else:
            cadd_phred = ''
            cadd_priPhCons = ''
            cadd_GerpRS = ''

        #current_cadd_str = cadd_phred + "\t" + v) for v in record.INFO['AF'])
        return cadd_phred, cadd_priPhCons, cadd_GerpRS



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

    outputfile.write("chr\tpos\tref\talt\tannotation\tgene_name\tlof" \
            "\texon\taa_pos\tpoly/sift\tAF\tGMAF\t1kgEMAF\tESPEMAF\t" \
            "HETEUR\tHOMEUR\tCADD\tpriPhCons\tGerpRS\n")


    vcf_reader = vcf.Reader(open(args.vcf, 'r'))
    for record in vcf_reader:
        current_chr = record.CHROM
        current_pos = record.POS
        current_ref = record.REF
        current_alt = ','.join(str(v) for v in record.ALT)
        current_af = ','.join(str(v) for v in record.INFO['AF'])
        current_het_nfe = ','.join(str(v) for v in record.INFO['Het_NFE'])
        current_hom_nfe = ','.join(str(v) for v in record.INFO['Hom_NFE'])

        # VEP
        if "CSQ" in record.INFO:
            csq = record.INFO['CSQ'][0].split('|')
            current_feature = csq[2]
            current_feature_type = csq[3]
            current_consequence = csq[4]
            current_sift = csq[24].split("(")[0]
            current_polyphen = csq[25].split("(")[0]
            current_eur_maf = csq[34]
            current_ea_maf = csq[37]
            current_LOF = csq[48]
            current_gmaf = csq[31]
        else:
            current_feature = ''
            current_feature_type = ''
            current_consequence = ''
            current_sift = ''
            current_polyphen = ''
            current_eur_maf = ''
            current_ea_maf = ''
            current_LOF = ''
            current_gmaf = ''


        # SnpEff
        ann = record.INFO['ANN'][0].split('|')
        annotation = ann[1]
        # annotation_impact = ann[2]
        current_gene = ann[3]
        current_exon = ann[8]
        current_aa_pos = ann[13]

        #CADD
        # cadd_phred, cadd_priPhCons, cadd_GerpRS
        (cadd_snp_phred, cadd_snp_priPhCons, cadd_snp_GerpRS) = getcadd(cadd_tbx, current_chr, current_pos)
        (cadd_indel_phred, cadd_indel_priPhCons, cadd_indel_GerpRS) = getcadd(cadd_indel_tbx, current_chr, current_pos)
        #for row in cadd_tbx.fetch(current_chr, current_pos-1, current_pos):
        #    row_info = row.split("\t")
        #    cadd_ref = row_info[2]
        #    cadd_alt = row_info[4]
        #    if(cadd_ref == current_ref and cadd_alt == current_alt):
        #        cadd_phred = row_info[115]
        #        cadd_priPhCons = row_info[18]
        #        cadd_GerpRS = row_info[26]
        #    else:
        #        cadd_phred = ''
        #        cadd_priPhCons = ''
        #        cadd_GerpRS = ''


        if "damaging" in current_polyphen or "deleterious" in current_sift:
            current_polysift = "del"
        else:
            current_polysift = ''


        out_str = [ "chr"+current_chr, str(current_pos), current_ref, current_alt,
                annotation, current_gene, current_LOF, current_exon,
                current_aa_pos, current_polysift, current_af, current_gmaf,
                current_eur_maf, current_ea_maf, current_het_nfe, current_hom_nfe,
                #cadd_snp, cadd_indel ]
                cadd_snp_phred, cadd_snp_priPhCons, cadd_snp_GerpRS ]
        out_str = [x or '.' for x in out_str]
        outputfile.write("\t".join(out_str))
        outputfile.write("\n")

    outputfile.close()



if __name__ == "__main__":
    main(sys.argv)

