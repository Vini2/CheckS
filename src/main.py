#!/usr/bin/env python3

import glob
import subprocess
import os
import re
import time
import argparse
import logging

from collections import defaultdict

from checks_utils import general_utils
from checks_utils import kmer_utils
from checks_utils import gene_utils
from checks_utils import plot_utils

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
reads_file_1 = args["reads1"]
reads_file_2 = args["reads2"]
query_ref_file = args["query_ref"]
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
logger.info("Forward reads file: "+reads_file_1)
logger.info("Reverse reads file: "+reads_file_2)
logger.info("Query reference genome file: "+query_ref_file)
logger.info("Folder with other reference genomes: "+other_refs)
logger.info("Extension of reference genome files: "+ref_ext)
logger.info("Final output folder: "+output_path)
logger.info("k value: "+str(k_val))
logger.info("Similarity threshold: "+str(similarity))
logger.info("Number of threads: "+str(nthreads))

logger.info("CheckS started")

start_time = time.time()


logger.info("Preparing files for the pipeline")

# Get Ref names

ref_files = glob.glob(other_refs+'*.'+ref_ext)
ref_files.sort()
ref_files += [query_ref_file]

query_ref_name = query_ref_file.split("/")[-1][:-(len(ref_ext)+1)]

contigs_file_ext = contigs_file.split(".")[-1]
contigs_file_name = contigs_file.split("/")[-1][:-(len(contigs_file_ext)+1)]

reads_file_ext = reads_file_1.split(".")[-1]
reads_file_1_name = reads_file_1.split("/")[-1][:-(len(reads_file_ext)+1)]
reads_file_2_name = reads_file_2.split("/")[-1][:-(len(reads_file_ext)+1)]

ref_names = general_utils.get_ref_names(ref_files, ref_ext)


# Make output folders for each ref and prepare files in each
   
output_folders = general_utils.make_output_folders(output_path,ref_names)
query_ref_folder = output_path+query_ref_name
       
general_utils.prep_files(ref_files, output_folders)


# K-mer analysis
logger.info("Running k-mer analysis")

## Run DSK
kmer_utils.run_dsk(output_folders, k_val, nthreads, output_path)


## Analyse k-mers for each ref

k_mer_res = {}

for i in range(len(output_folders)):
    
    logger.info("Processing files in folder: "+output_folders[i])

#     Get ref k-mers
    ref_kmers = {}

    with open(output_folders[i]+"/ref.txt", "r") as myfile:

        for line in myfile.readlines():
            strings = line.strip().split(" ")
            kmer = strings[0]
            kmer_count = int(strings[1])

            ref_kmers[kmer] = kmer_count
            
    logger.debug("Total number of k-mers of "+ref_names[i]+": "+str(len(ref_kmers)))
            
#     Get k-mers of other references without ref
    other_kmers = []

    with open(output_folders[i]+"/refs_combined.txt", "r") as myfile:

        for line in myfile.readlines():
            strings = line.strip().split(" ")
            kmer = strings[0]
            other_kmers.append(kmer)

    other_kmers = set(other_kmers)
    
    
#     Get unique k-mers of ref
    unique_kmers = set(ref_kmers.keys()).difference(other_kmers)
    
    logger.debug("Unique k-mers of "+ref_names[i]+": "+str(len(unique_kmers)))
    k_mer_res[ref_names[i]] = [len(unique_kmers)]
    
    
#     Get k-mers common in ref and reads
    common_read_kmers = {}

    with open(output_path+"/reads.txt", "r") as myfile:

        for line in myfile.readlines():

            strings = line.strip().split(" ")
            kmer = strings[0]
            kmer_count = int(strings[1])

            if kmer in unique_kmers:
                common_read_kmers[kmer] = kmer_count

    logger.debug("Unique k-mers of "+ref_names[i]+" found in reads: "+str(len(common_read_kmers)))
    k_mer_res[ref_names[i]].append(len(common_read_kmers))
    
    
#     Get k-mers common in ref and contigs
    common_contig_kmers = {}

    with open(output_path+"contigs.txt", "r") as myfile:

        for line in myfile.readlines():

            strings = line.strip().split(" ")
            kmer = strings[0]
            kmer_count = int(strings[1])

            if kmer in unique_kmers:
                common_contig_kmers[kmer] = kmer_count

    logger.debug("Unique k-mers of "+ref_names[i]+" found in contigs: "+str(len(common_contig_kmers)))
    k_mer_res[ref_names[i]].append(len(common_contig_kmers))


## Plot k-mer analysis results

kmers_in_contigs = []
kmers_in_reads = []

for ref in k_mer_res:
    kmers_in_reads.append(k_mer_res[ref][1]/k_mer_res[ref][0]*100)
    kmers_in_contigs.append(k_mer_res[ref][2]/k_mer_res[ref][0]*100)

plot_utils.plot_kmer_res(0.5, 20, 15, 'png', 300, kmers_in_reads, kmers_in_contigs, k_mer_res.keys(), output_path, prefix)


# Coverage analysis
logger.info("Running coverage analysis")

## Run CoverM
if not os.path.isfile(output_path+"/coverage_output.tsv"):
    subprocess.run("coverm genome --coupled "+reads_file_1+" "+reads_file_2+" --genome-fasta-files "+" ".join(ref_files)+" -o "+output_path+"/coverage_output.tsv --threads "+str(n_threads), shell=True)
else:
    logger.debug("CoverM result exists.")

logger.info("CoverM results can be found in "+output_path+"/coverage_output.tsv")


# Assembly analysis
logger.info("Running assembly analysis")

## Run metaQUAST

if not os.path.isdir(output_path+"quast_results"):
    logger.info("Running metaQUAST")
    subprocess.run("metaquast -r "+",".join(ref_files)+" -o "+output_path+"quast_results -t "+str(8)+" "+contigs_file, shell=True)

logger.info("The metaQUAST results can be found in "+output_path+"quast_results")



# Gene analysis
logger.info("Running gene analysis")

## Run MetaGeneMark
gene_utils.run_metagenemark(output_folders, contigs_file, contigs_file_name, contigs_file_ext, output_path)


## Analyse genes for each ref

gene_res = {}

for i in range(len(output_folders)):
    
    logger.info("Processing files in folder "+output_folders[i])
    
    # Get ref gene sequences
    ref_gene_nt_seq = []

    gene_nt_seq = {}

    for seq in gene_utils.get_seqs(output_folders[i]+"/ref.fasta.lst"):
        if "_nt|" in seq:

            seqs = seq.split("\n")

            my_seq = ""

            for n in range(len(seqs)):
                if n > 0:
                    my_seq += seqs[n]

            strings = seq.split("|")

            gene_id = strings[0][1:]

            gene_nt_seq[gene_id] = seq

            ref_gene_nt_seq.append(my_seq)

    with open(output_folders[i]+"ref_genes.fna", "w") as ntfile:
        for gene in gene_nt_seq:
            ntfile.write(gene_nt_seq[gene])

    logger.debug("Total number of genes found in the reference: "+str(len(ref_gene_nt_seq)))
    
    
    # Get gene sequences of remaining references
    other_gene_nt_seq = []
    gene_nt_seq = {}

    for seq in gene_utils.get_seqs(output_folders[i]+"/refs_combined.fasta.lst"):
        if "_nt|" in seq:

            seqs = seq.split("\n")

            my_seq = ""

            for n in range(len(seqs)):
                if n > 0:
                    my_seq += seqs[n]

            strings = seq.split("|")

            gene_id = strings[0][1:]

            gene_nt_seq[gene_id] = seq

            other_gene_nt_seq.append(my_seq)

    with open(output_folders[i]+"refs_combined_genes.fna", "w") as ntfile:
        for gene in gene_nt_seq:
            ntfile.write(gene_nt_seq[gene])

    logger.debug("Total number of genes found in other reference: "+str(len(other_gene_nt_seq)))
    
    
    # Get common genes between ref and other refs
    
    gene_utils.align(output_folders[i]+"refs_combined_genes.fna", output_folders[i]+"ref_genes.fna", output_folders[i], "refs_mapped_genes.output", 0.9)
    
    commons_genes = []

    with open(output_folders[i]+"refs_mapped_genes.output", 'r') as f:
        for line in f.readlines():
            strings = line.strip().split()

            if strings[1] != "POOR_MAPPING":

                gene_id = strings[0].split("|")[0]

                commons_genes.append(gene_id)

    uniq_gene_count = len(ref_gene_nt_seq)-len(commons_genes)
    
    logger.debug("Number of unique reference genes: "+str(uniq_gene_count))
    gene_res[ref_names[i]] = [uniq_gene_count]
    
    
    gene_nt_seq = {}

    for seq in gene_utils.get_seqs(output_folders[i]+"/ref.fasta.lst"):
        if "_nt|" in seq:

            seqs = seq.split("\n")

            my_seq = ""

            for n in range(len(seqs)):
                if n > 0:
                    my_seq += seqs[n]

            strings = seq.split("|")

            gene_id = strings[0][1:]

            if gene_id not in commons_genes:
                gene_nt_seq[gene_id] = seq


    with open(output_folders[i]+"ref_unique_genes.fna", "w") as ntfile:
        for gene in gene_nt_seq:
            ntfile.write(gene_nt_seq[gene])

    
    # Get gene sequences of contigs
    contig_genes_file = output_path+contigs_file_name+"."+contigs_file_ext+".lst"
    
    contig_nt_seq = []

    gene_nt_seq = {}

    for seq in gene_utils.get_seqs(contig_genes_file):
        if "_nt|" in seq:

            seqs = seq.split("\n")

            my_seq = ""

            for n in range(len(seqs)):
                if n > 0:
                    my_seq += seqs[n]

            strings = seq.split("|")

            gene_id = strings[0][1:]

            gene_nt_seq[gene_id] = seq

            contig_nt_seq.append(my_seq)

    if not os.path.isfile(output_path+"contigs_genes.fna"):
        with open(output_path+"contigs_genes.fna", "w") as ntfile:
            for gene in gene_nt_seq:
                ntfile.write(gene_nt_seq[gene])
    
    
    # Get common genes between ref unique genes and contig genes
    gene_utils.align(output_path+"contigs_genes.fna", output_folders[i]+"ref_unique_genes.fna", output_folders[i], "refs_unique_contigs_mapped_genes.output", 0.9)
    
    commons_genes = []

    with open(output_folders[i]+"refs_unique_contigs_mapped_genes.output", 'r') as f:
        for line in f.readlines():
            strings = line.strip().split()

            if strings[1] != "POOR_MAPPING":

                gene_id = strings[0].split("|")[0]

                commons_genes.append(gene_id)

    logger.debug("Unique ref genes mapped to contigs: "+str(len(commons_genes)))
    gene_res[ref_names[i]].append(len(commons_genes))


## Plot gene analysis results

genes_in_contigs = []

for ref in gene_res:
    genes_in_contigs.append(gene_res[ref][1]/gene_res[ref][0]*100)

plot_utils.plot_gene_res(0.5, 20, 15, 'png', 300, genes_in_contigs, gene_res.keys(), output_path, prefix)


# Get elapsed time

## Determine elapsed time
elapsed_time = time.time() - start_time

## Print elapsed time for the process
logger.info("Elapsed time: "+str(elapsed_time)+" seconds")



# Exit program

logger.info("Thank you for using CheckS!")