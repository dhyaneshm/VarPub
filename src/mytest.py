import pybedtools

a = pybedtools.example_bedtool('a.bed')
b = pybedtools.example_bedtool('b.bed')

print "cat a.bed\n" + str(a)
print "cat b.bed\n" + str(b)

print a.intersect(b)

