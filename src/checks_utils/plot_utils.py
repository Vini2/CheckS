#!/usr/bin/env python3

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker
import logging

plt.rcParams.update({'font.size': 24})

# create logger
logger = logging.getLogger('CheckS 0.1')


def plot_kmer_res(width, size1, size2, figformat, figdpi, kmers_in_reads, kmers_in_contigs, x_labels, output_path, prefix):

    mytitle = "Unique reference k-mers found in reads and contigs"

    ind = np.arange(len(kmers_in_reads))  # the x locations for the groups

    fig, ax = plt.subplots(figsize=(size1, size2))
    rects1 = ax.bar(ind-width/2, kmers_in_reads, width,
                    yerr=(0), label='Common k-mers in reads')
    rects2 = ax.bar(ind+width/2, kmers_in_contigs, width,
                    yerr=(0), label='Common k-mers in contigs')

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Precentage found (%)', labelpad=0)
    ax.set_xlabel("Reference genomes")
    ax.set_ylim((0, 105))
    ax.set_title(mytitle)
    ax.set_xticks(ind)
    ax.set_xticklabels(x_labels, rotation=90)
    ax.legend(loc='upper right', framealpha=1, prop={'size': 20})
    ax.set_axisbelow(True)
    ax.yaxis.grid(True)

    fig.tight_layout()

    plt.savefig(output_path+prefix+'kmer_res.'+figformat,
                format=figformat, dpi=figdpi, bbox_inches='tight')
    logger.info("Plot of k-mer analysis results can be found in " +
                output_path+prefix+'kmer_res.'+figformat)


def plot_gene_res(width, size1, size2, figformat, figdpi, genes_in_contigs, x_labels, output_path, prefix):

    mytitle = "Unique reference genes found in contigs"

    ind = np.arange(len(genes_in_contigs))  # the x locations for the groups

    fig, ax = plt.subplots(figsize=(size1, size2))
    rects1 = ax.bar(ind, genes_in_contigs, width, yerr=(0),
                    label='Common genes in contigs')

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Precentage found (%)', labelpad=0)
    ax.set_xlabel("Reference genomes")
    ax.set_ylim((0, 105))
    ax.set_title(mytitle)
    ax.set_xticks(ind)
    ax.set_xticklabels(x_labels, rotation=90)
    ax.set_axisbelow(True)
    ax.yaxis.grid(True)

    fig.tight_layout()

    plt.savefig(output_path+prefix+'gene_res.'+figformat,
                format=figformat, dpi=figdpi, bbox_inches='tight')
    logger.info("Plot of gene analysis results can be found in " +
                output_path+prefix+'gene_res.'+figformat)
