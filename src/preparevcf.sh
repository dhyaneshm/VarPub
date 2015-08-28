#!/bin/bash
# Author: Khalid Mahmood
# Prepare a GATK-best practices based VCF file
# for downstream analysis.

# Some variables

VCF=kr.vcf
VCFNAME=${VCF%.*}
NORMVCF=$VCFNAME.norm.vcf
VEPVCF=$VCFNAME.norm.vep.vcf
SNPEFFVCF=$VCFNAME.norm.vep.snpeff.vcf

SNPS=$VCFNAME.snps.vcf
INDELS=$VCFNAME.indels.vcf

GATKSNPS=$VCFNAME.snps.gatk.vcf
GATKSNPSPASS=$VCFNAME.snps.gatk.pass.vcf
GATKINDELS=$VCFNAME.indels.gatk.vcf
GATKINDELSPASS=$VCFNAME.indels.gatk.pass.vcf

REF=genome.fa
VEPPATH=/vlsci/VR0002/kmahmood/Programs/vep/77/ensembl-tools-release-77/scripts/variant_effect_predictor
SNPEFFJAR=/usr/local/snpeff/4.1g/snpEff.jar
SNPEFFCONF=/vlsci/VR0002/kmahmood/reference/snpeff/4.1/snpEff.config
GATK=/usr/local/gatk/3.2-2/GenomeAnalysisTK.jar

AF=0.11
EXACAF=0.01

# Decompose ... etc.
zless $VCF \
    | sed 's/ID=AD,Number=./ID=AD,Number=R/' \
    | vt decompose -s - \
    | vt normalize -r $REF - > $NORMVCF
wait
#    | bgzip -c > $NORMVCF
# tabix -p vcf $NORMVCF

# Run VEP annotaion - includes ExAC_AN and ExAC_AC numbers for AF calculation
perl $VEPPATH/variant_effect_predictor.pl --cache -i $NORMVCF --format vcf -o $VEPVCF --force_overwrite --vcf --fork 4 --everything --offline --coding_only --no_intergenic --plugin LoF,human_ancestor_fa:~/.vep/homo_sapiens/77_GRCh37/human_ancestor.fa.gz -custom ~/reference/ExAC0.3/ExAC.r0.3.sites.vep.vcf.gz,ExAC,vcf,exact,0,AF,AC,AC_AFR,AC_AMR,AC_Adj,AC_EAS,AC_FIN,AC_Het,AC_Hom,AC_NFE,AC_OTH,AC_SAS,AF,AN,AN_AFR,AN_AMR,AN_Adj,AN_EAS,AN_FIN,AN_NFE,AN_OTH,AN_SAS
wait

# RUN SnpEff ANNOTATION
java -jar $SNPEFFJAR eff -c $SNPEFFCONF -canon hg19 $VEPVCF > $SNPEFFVCF
wait

# SPLIT IN TO SNPs AND INDELs
vt view -h $SNPEFFVCF -f "VTYPE==SNP&&N_ALLELE==2&&INFO.AF<$AF" > $SNPS
vt view -h $SNPEFFVCF -f "VTYPE==SNP&&N_ALLELE==2&&INFO.AF<$AF" > $INDELS
wait

# RUN GATK SNPs + INDELS
zless $SNPS | sed 's/ID=AD,Number=R/ID=AD,Number=./' > $SNPS.temp
wait
mv $SNPS.temp $SNPS
java -Xmx4g -jar $GATK -R $REF -T VariantFiltration -V $SNPS -o $GATKSNPS --clusterWindowSize 10 --filterExpression "QUAL < 30.0 " --filterName "VeryLowQual" --filterExpression "QD < 2.0 " --filterName "LowQD" --filterExpression "DP < 10 " --filterName "LowCoverage" --filterExpression "MQ < 40 " --filterName "LowMappingQual"  --filterExpression "FS > 60.0 " --filterName "StrandBias"
wait
zless $INDELS | sed 's/ID=AD,Number=R/ID=AD,Number=./' > $INDELS.temp
mv $INDELS.temp $INDELS
wait
java -Xmx4g -jar $GATK -R $REF -T VariantFiltration -V $INDELS -o $GATKINDELS --clusterWindowSize 10 --filterExpression "QUAL < 30.0 " --filterName "VeryLowQual" --filterExpression "QD < 2.0 " --filterName "LowQD" --filterExpression "DP < 10 " --filterName "LowCoverage" --filterExpression "MQ < 40 " --filterName "LowMappingQual" --filterExpression "FS > 200.0 " --filterName "StrandBias"
wait

# GATK PASS ONLY
zless $GATKSNPS | sed 's/ID=AD,Number=./ID=AD,Number=R/' > $GATKSNPS.temp
mv $GATKSNPS.temp $GATKSNPS
zless $GATKINDELS | sed 's/ID=AD,Number=./ID=AD,Number=R/' > $GATKINDELS.temp
mv $GATKINDELS.temp $GATKINDELS
vt view -h $GATKSNPS -f "PASS" > $GATKSNPSPASS
vt view -h $GATKINDELS -f "PASS" > $GATKINDELSPASS
zless $GATKSNPSPASS | sed 's/ID=AD,Number=R/ID=AD,Number=./' > $GATKSNPS.temp
mv $GATKSNPS.temp $GATKSNPSPASS
zless $GATKINDELSPASS | sed 's/ID=AD,Number=R/ID=AD,Number=./' > $GATKINDELS.temp
mv $GATKINDELS.temp $GATKINDELSPASS

