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

def getcadd(cadd_tbx, current_chr, current_pos, current_ref, current_alt):
    data = cadd_tbx.fetch(current_chr, current_pos-1, current_pos)
    cadd_phred, cadd_priPhCons, cadd_GerpRS = '','',''
    cadd_polysift, cadd_test1, cadd_test2 = '','',''
    if data is not None:
        for row in data:
            row_info = row.split("\t")
            cadd_ref = row_info[2]
            cadd_alt = row_info[4]
            if(cadd_ref == current_ref and cadd_alt == current_alt):
                cadd_phred = row_info[115]
                cadd_priPhCons = row_info[18]
                cadd_GerpRS = row_info[26]
                if "damaging" in row_info[110] or "deleterious" in row_info[112]:
                    cadd_polysift = "del"
                break
    #else:
    #    cadd_phred = ''
    #    cadd_priPhCons = ''
    #    cadd_GerpRS = ''

    #current_cadd_str = cadd_phred + "\t" + v) for v in record.INFO['AF'])
    return cadd_phred, cadd_priPhCons, cadd_GerpRS, \
            cadd_polysift

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
        #current_alt_array = current_alt.split(",")
        current_af = ','.join(str(v) for v in record.INFO['AF'])
        current_het_nfe = ','.join(str(v) for v in record.INFO['Het_NFE'])
        current_hom_nfe = ','.join(str(v) for v in record.INFO['Hom_NFE'])
        # CHECK INDEL AND MNP
        indel = True if ( len(current_ref) > 1 or len(current_alt) > 1 and ("," not in current_alt ) ) else False
        # mnp = map(labmda x, len(record.ALT)
        mnp = True if len(record.ALT) > 1 else False

        #for current_alt in record.ALT:

        # VEP
        if "CSQ" in record.INFO:
            csq = record.INFO['CSQ'][0].split('|')
            current_feature, current_feature_type = csq[2], csq[3]
            current_consequence = csq[4]
            current_sift = csq[24].split("(")[0]
            current_polyphen = csq[25].split("(")[0]
            current_gmaf, current_eur_maf, current_ea_maf =  csq[31], csq[34], csq[37]
            current_LOF = csq[48]
        else:
            current_feature, current_feature_type, current_consequence = '','',''
            current_sift, current_polyphen, current_eur_maf = '','',''
            current_ea_maf, current_LOF, current_gmaf = '','',''


        # SnpEff
        ann = record.INFO['ANN'][0].split('|')
        annotation = ann[1]
        #   GENE INFORMATION
        current_gene, current_exon, current_aa_pos = ann[3], ann[8], ann[13]

        #CADD SNP
        cadd_phred_temp = ''
        cadd_phred = ''
        indel = ''
        mnp_cadds = []
        if not indel:
            for alt in record.ALT:
                (cadd_phred_temp, cadd_snp_priPhCons, cadd_snp_GerpRS, cadd_polysift) = \
                        getcadd(cadd_tbx, current_chr, current_pos, current_ref, alt)
                mnp_cadds.append(str(alt) + ":" + cadd_phred_temp)
                cadd_phred = ",".join(mnp_cadds)
                indel = "F"
        else:
            #CADD INDEL
            (cadd_phred, cadd_indel_priPhCons, cadd_indel_GerpRS, cadd_polysift) = \
                    getcadd(cadd_indel_tbx, current_chr, current_pos, current_ref, current_alt)
            indel = "T"

        #if cadd_snp_phred:
        #    single_cadd_score = cadd_snp_phred
        #else:
        #    single_cadd_score = cadd_indel_phred


        out_str = [ "chr"+current_chr, str(current_pos), current_ref, current_alt,
                annotation, current_gene, current_LOF, current_exon,
                current_aa_pos, cadd_polysift, current_af, current_gmaf,
                current_eur_maf, current_ea_maf, current_het_nfe, current_hom_nfe,
                cadd_phred, cadd_snp_priPhCons, cadd_snp_GerpRS, indel ]
        out_str = [x or '.' for x in out_str]
        outputfile.write("\t".join(out_str))
        outputfile.write("\n")
        #cadd_polysift = cadd

    outputfile.close()


if __name__ == "__main__":
    main(sys.argv)

