# -*- coding: utf-8 -*-
"""
Date: Mon Aug 21 12:51:16 2023
File Name: Run_Pipeline.py
Author: Mirjam Karlsson-Müller

Description:
    Framework for the pipeline. Reads in inputs, allocates resources and runs 
    different parts of the pipeline.
    
List of Functions:
    
Procedure: 
    1.
    2.
    3.
    
Useage:
    
    
Possible Bugs:
"""

#%% 0.0 Imports

import argparse
import re
import time
import multiprocessing as mp

#%% 0.1 Argparse

parser = argparse.ArgumentParser(prog='Score Alternative Splicing',
                                 usage='%(prog)s -s SAMPLE-FOLDER -o OUTPUT-FOLDER \
                                     -g GENCODE-FILE -r REFSEQ-FILE \
                                         [-c] "chrX:XXXXXX-XXXXXX" -as AS-TYPE \
                                             -is INSERT-SIZE',
                                 description="""Creates
                                 a table with the PSI scores supporting said event
                                 per sample in sample folder.""")

parser.add_argument('--samples', '-s', required=True,
                    help='file containing the names of all samples that the pipeline is to be run for.')
parser.add_argument('--genes', '-g', required=True,
                    help="File containing all HGNC gene names, the pipeline should run for.")

args = parser.parse_args()


#%% 0.2 Functions





#%% 0.3 Time Start

start_time=time.time()

#%% 1. Read in Samples and Genes


#%% 2. Initialize Pipeline for every Gene, Allocate Pararell resources!









#%% Time End

print("Run time: {:.2f} seconds.".format(time.time()-start_time))  