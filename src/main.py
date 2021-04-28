#!/usr/bin/env python3

import sys
import csv
import time
import argparse
import re
import heapq
import os
import math
import networkx as nx
import logging
import numpy as np
import operator
import gc
import subprocess

from multiprocessing import Pool
from Bio import SeqIO
from igraph import *
from collections import defaultdict
from tqdm import tqdm
from functools import partial

# Set paramters
#---------------------------------------------------

# MAX_WEIGHT = sys.float_info.max
# MG_LENGTH_THRESHOLD = 0.3
# SEED_MG_THRESHOLD = 0.1
# DEPTH = 1


# Setup argument parser
#---------------------------------------------------

parser = argparse.ArgumentParser(description="""CheckS is a pipeline to check the presence of a species
                                            in a metagenomics sample, provided its reference genome and
                                            the reference genomes of the other known species.""")

parser.add_argument("--contigs", 
                    required=True,
                    type=str,
                    help="path to the contigs file")

parser.add_argument("--reads1", 
                    required=True,
                    type=str,
                    help="path to the forward reads file")

parser.add_argument("--reads2", 
                    required=True,
                    type=str,
                    help="path to the reverse reads file")

parser.add_argument("--query_ref", 
                    required=True,
                    type=str,
                    help="path to the query reference genome")

parser.add_argument("--other_refs", 
                    required=True,
                    type=str,
                    help="path to the folder with other reference genomes")

parser.add_argument("--ref_ext", 
                    required=True, 
                    type=str, 
                    default="fasta", 
                    help="extension of the reference genome files. [default: fasta]")

parser.add_argument("--output", 
                    required=True,
                    type=str,
                    help="path to the output folder")

parser.add_argument("--k", 
                    required=False, 
                    type=int, 
                    default=25, 
                    help="k value to run DSK k-mer counting. [default: 25]")

parser.add_argument("--similarity", 
                    required=False, 
                    type=float, 
                    default=0.9, 
                    help="similarity threshold for mapping. [default: 0.9]")

parser.add_argument("--prefix", 
                    required=False,
                    type=str,
                    default='',
                    help="prefix for the output file")

parser.add_argument("--nthreads", 
                    required=False, 
                    type=int, 
                    default=8, 
                    help="number of threads/cores to use. [default: 8]")

args = vars(parser.parse_args())

contigs_file = args["contigs"]
reads1 = args["reads1"]
reads2 = args["reads2"]
query_ref = args["query_ref"]
other_refs = args["other_refs"]
ref_ext = args["ref_ext"]
output_path = args["output"]
k_val = args["k"]
similarity = args["similarity"]
prefix = args["prefix"]
nthreads = args["nthreads"]


# Setup logger
#-----------------------
logger = logging.getLogger('CheckS 0.1')
logger.setLevel(logging.DEBUG)
logging.captureWarnings(True)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
consoleHeader = logging.StreamHandler()
consoleHeader.setFormatter(formatter)
consoleHeader.setLevel(logging.INFO)
logger.addHandler(consoleHeader)

# Setup output path for log file
#---------------------------------------------------

fileHandler = logging.FileHandler(output_path+"/"+prefix+"checks.log")
fileHandler.setLevel(logging.DEBUG)
fileHandler.setFormatter(formatter)
logger.addHandler(fileHandler)

logger.info("Welcome to CheckS: Checking the presence of a species in a metagenomics sample.")

logger.info("Input arguments:")
logger.info("Contigs file: "+contigs_file)
logger.info("Forward reads file: "+reads1)
logger.info("Reverse reads file: "+reads2)
logger.info("Query reference genome file: "+query_ref)
logger.info("Folder with other reference genomes: "+other_refs)
logger.info("Extension of reference genome files: "+ref_ext)
logger.info("Final output folder: "+output_path)
logger.info("k value: "+str(k_val))
logger.info("Similarity threshold: "+str(similarity))
logger.info("Number of threads: "+str(nthreads))

logger.info("CheckS started")

start_time = time.time()