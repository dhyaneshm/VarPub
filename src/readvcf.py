import pybedtools
import os

cpath = os.path.dirname(os.path.abspath(__file__))

kgbed = pybedtools.BedTool('/Users/khalidm/Development/repositories/VarPub/data/1kg.phase1.snp.bed.gz')
avcf = pybedtools.BedTool('/Users/khalidm/Development/repositories/VarPub/data/BRCA2.snpeff.vcf')

# print "cat a.bed\n" + str(avcf)
# print "cat b.bed\n" + str(bfile)

# print intersect
# print avcf.intersect(kgbed)
# print length
# print len(a)

print avcf[0][7] + " -> " + avcf[0][8]

