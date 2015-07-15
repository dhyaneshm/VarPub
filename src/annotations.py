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

def getTabixVal(input_tbx, current_chr, current_pos, current_ref, current_alt):
    #current_chr = current_chr.translate(None, 'chr')
    data = input_tbx.fetch(current_chr, current_pos-1, current_pos)
    value = '.'

    if data is not None:
        for row in data:
            row_info = row.split("\t")
            value = row_info[3]
    else:
        value = '.'

    return value


def getTabixBool(input_tbx, current_chr, current_pos, current_ref, current_alt):
    #current_chr = current_chr.translate(None, 'chr')
    data = input_tbx.fetch(current_chr, current_pos-1, current_pos)
    value = '.'

    if data is not None:
        return True
    else:
        return False


