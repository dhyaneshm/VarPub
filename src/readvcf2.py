'''
A tool to annotate and print variants in tabular format.
Author: Khalid Mahmood (khalid.mahmood@unimelb.edu.au).
Copyright: 2015
'''

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
    current_chr = current_chr.translate(None, 'chr')
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
    else:
        cadd_phred = '.'

    return cadd_phred, cadd_priPhCons, cadd_GerpRS, \
            cadd_polysift

def getfathmm(fathmm_tbx, current_chr, current_pos, current_ref, current_alt):
    current_chr = current_chr.translate(None, 'chr')
    data = fathmm_tbx.fetch(current_chr, current_pos-1, current_pos)
    fathmm_score = ''
    if data is not None:
        for row in data:
            row_info = row.split("\t")
            fathmm_ref = row_info[3]
            fathmm_alt = row_info[4]
            if(fathmm_ref == current_ref and fathmm_alt == current_alt):
                fathmm_score = row_info[7]
                break
    # else:
    #    fathmm_score = ''

    return fathmm_score

# return allele frequency given the allele count and assuming allele number = (total allele number/2)
def getAF(ac, an):
    if(float(an)>0):
        af_temp = ac / an
        newlist = round(af_temp, 5)
    else:
        newlist = 'NA'
    return str(newlist)

# return index of the current alt allele from exac multiallelic data
def getexacallele(exac_tbx, current_chr, current_pos, current_ref, current_alt):
    current_chr = current_chr.translate(None, 'chr')
    data = exac_tbx.fetch(current_chr, current_pos-1, current_pos)
    index = -2
    exac_row = next(data, None)
    if exac_row:
        exac_alt_temp = exac_row.split("\t")[4]
        exac_alt_row = exac_alt_temp.split(",")
        indexlist = [i for i, s in enumerate(exac_alt_row) if current_alt is s]
        index = indexlist[-1] if len(indexlist)==1 else -2
        #print "\t\tT\t" + str(index)
    else:
        index = -2
        #print "\t\tF\t" + str(index)
    return index


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
        current_chr = record.CHROM
        current_id = record.ID
        current_pos = record.POS
        current_ref = record.REF
        current_alt = ','.join(str(v) for v in record.ALT)
        #current_alt_array = current_alt.split(","
        current_af = ','.join(str(v) for v in record.INFO['AF'])
        current_het_nfe = ''
        current_hom_nfe = ''

        # check if the variant is in ExAC annotated
        if any("ExAC" in s for s in record.INFO):
            #print current_chr + "\t" + current_id + "\t" + current_ref + ":" + current_alt + str(record.INFO['ExAC_AN_Adj']) + "\t" + str(record.INFO['ExAC_AN_Adj'])
            current_exac_index = getexacallele(exac_tbx, current_chr, current_pos, current_ref, current_alt)
            if(current_exac_index>-2):
                current_het_nfe = ','.join(str(v) for v in record.INFO['ExAC_AC_Het'])
                current_hom_nfe = ','.join(str(v) for v in record.INFO['ExAC_AC_Hom'])
                current_exac_af = getAF(float(record.INFO['ExAC_AC_Adj'][current_exac_index]),float(record.INFO['ExAC_AN_Adj'][-1])) # Total adjusted
                current_exac_eas = getAF(float(record.INFO['ExAC_AC_EAS'][current_exac_index]),float(record.INFO['ExAC_AN_EAS'][-1])) # East Asians
                current_exac_nfe = getAF(float(record.INFO['ExAC_AC_NFE'][current_exac_index]),float(record.INFO['ExAC_AN_NFE'][-1])) # NonFin Eur
                current_exac_fin = getAF(float(record.INFO['ExAC_AC_FIN'][current_exac_index]),float(record.INFO['ExAC_AN_FIN'][-1])) # Fin Eur
                current_exac_sas = getAF(float(record.INFO['ExAC_AC_SAS'][current_exac_index]),float(record.INFO['ExAC_AN_SAS'][-1])) # South Asian
                current_exac_afr = getAF(float(record.INFO['ExAC_AC_AFR'][current_exac_index]),float(record.INFO['ExAC_AN_AFR'][-1])) # African
                current_exac_amr = getAF(float(record.INFO['ExAC_AC_AMR'][current_exac_index]),float(record.INFO['ExAC_AN_AMR'][-1])) # Latino
                current_exac_oth = getAF(float(record.INFO['ExAC_AC_OTH'][current_exac_index]),float(record.INFO['ExAC_AN_OTH'][-1])) # Other
        else:
            current_exac_af,current_exac_eas,current_exac_nfe = 0,0,0
            current_exac_fin,current_exac_sas,current_exac_afr = 0,0,0
            current_exac_amr,current_exac_oth = 0,0

        # CHECK INDEL AND MNP
        #print current_ref + ":" + current_alt
        indel = True if ((len(current_ref) > 1 or len(current_alt) > 1) and \
                ("," not in current_ref and "," not in current_alt)) else False
        # mnp = map(labmda x, len(record.ALT)
        mnp = True if len(record.ALT) > 1 else False

        #for current_alt in record.ALT:

        # VEP
        current_sift, current_polyphen, current_consequence, current_LOF = '','','',''
        current_gmaf, current_eur_maf, current_ea_maf = '','',''
        current_feature, current_feature_type = '',''
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
        indel_str= ''
        mnp_cadds = []
        cadd_scores = []
        fathmm_score = ''
        for alt in record.ALT:
            if(len(current_ref) == 1 and len(alt) == 1):
                (cadd_phred_temp, cadd_priPhCons, cadd_GerpRS, cadd_polysift) = \
                        getcadd(cadd_tbx, current_chr, current_pos, current_ref, alt)
                mnp_cadds.append(str(alt) + ":" + cadd_phred_temp)
                cadd_scores.append(cadd_phred_temp)
                # GET FATHMM SCORE
                fathmm_score = getfathmm(fathmm_tbx, current_chr, current_pos, current_ref, alt)
            else: # IF VAR IS AN INDEL
                (cadd_phred_temp, cadd_priPhCons, cadd_GerpRS, cadd_polysift) = \
                        getcadd(cadd_indel_tbx, current_chr, current_pos, current_ref, alt)
                mnp_cadds.append(str(alt) + ":" + cadd_phred_temp)
                cadd_scores.append(cadd_phred_temp)
        cadd_phred = ",".join(mnp_cadds)
        # indel_str = "."

        out_str = [ current_chr, str(current_pos), str(current_id), current_ref, current_alt,
                annotation, current_gene, current_LOF, current_exon,
                current_aa_pos, cadd_polysift, current_af, current_gmaf,
                current_eur_maf, current_ea_maf,
                #current_het_nfe, current_hom_nfe,
                current_exac_af, current_exac_eas, current_exac_nfe, current_exac_fin, current_exac_sas,
                current_exac_afr, current_exac_amr, current_exac_oth,
                cadd_phred, str(max(cadd_scores)), cadd_priPhCons, cadd_GerpRS,
                fathmm_score ]
        out_str = [x or '.' for x in out_str]
        outputfile.write("\t".join(out_str))
        outputfile.write("\n")

    outputfile.close()


if __name__ == "__main__":
    main(sys.argv)

