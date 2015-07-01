#!/usr/bin/python

import sys
import os
import pybedtools

def findlist(input_list, a):
    result = []
    index = -2
    for i, x in enumerate(input_list):
        if x == a:
            index = i
            break
    return index

