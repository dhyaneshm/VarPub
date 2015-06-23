#!/usr/bin/python

import sys
import os
import argparse
import getopt
import pybedtools


class Config(object):
    def __init__(self, current_record, config_filename):
        # Try to open and parse the YAML formatted config file
        with open(config_filename) as config_file:
            try:
                config = yaml.load(config_file)
            except yaml.YAMLError, exc:
                print("Error in configuration file:", exc)
                raise exc
        self.config = config
        self.config_filename = config_filename





