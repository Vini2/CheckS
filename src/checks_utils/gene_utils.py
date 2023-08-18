#!/usr/bin/env python3

import subprocess
import os
import logging

from collections import defaultdict

# create logger
logger = logging.getLogger('CheckS 0.1.0')


def run_metagenemark(output_folders, contigs_file, contigs_file_name, contigs_file_ext, output_path):

    mgm_path = os.environ['MGM_PATH'] + "/"

    for output_folder in output_folders:

        logger.debug("Running MetaGeneMark in "+output_folder)

        # Run MetaGeneMark on ref.fasta
        subprocess.run(mgm_path + "gmhmmp -a -d -m " + mgm_path +
                       "MetaGeneMark_v1.mod "+output_folder+"ref.fasta", shell=True)

        # Run MetaGeneMark on refs_combined.fasta
        subprocess.run(mgm_path + "gmhmmp -a -d -m " + mgm_path +
                       "MetaGeneMark_v1.mod "+output_folder+"refs_combined.fasta", shell=True)

    # Run MetaGeneMark on contigs
    subprocess.run(mgm_path + "gmhmmp -a -d -m " + mgm_path +
                   "MetaGeneMark_v1.mod "+contigs_file, shell=True)
    subprocess.run("mv "+contigs_file+".lst "+output_path +
                   contigs_file_name+"."+contigs_file_ext+".lst", shell=True)


def get_seqs(path):
    active = ""

    for line in open(path):
        if line[0] == ">":
            active += line
        elif len(active) > 0 and len(line.strip()) != 0:
            active += line
        elif len(line.strip()) == 0 and len(active) > 0:
            yield active
            active = ""
    if len(active) > 0:
        yield active


def align(target, query, output_folder, output_name, threshold):

    if not os.path.isfile(output_folder+output_name):

        subprocess.run("minimap2 "+target+" "+query+" > " +
                       output_folder+output_name+".paf", shell=True)

        contig_ref = defaultdict(list)
        contig_ref_aln_length = defaultdict(list)
        contig_length = {}

        for line in open(output_folder+output_name+".paf"):
            data = line.strip().split('\t')

#     1	string	Query sequence name
#     2	int	Query sequence length
#     3	int	Query start (0-based; BED-like; closed)
#     4	int	Query end (0-based; BED-like; open)
#     5	char	Relative strand: "+" or "-"
#     6	string	Target sequence name
#     7	int	Target sequence length
#     8	int	Target start on original strand (0-based)
#     9	int	Target end on original strand (0-based)
#     10	int	Number of residue matches
#     11	int	Alignment block length
#     12	int	Mapping quality (0-255; 255 for missing)

            qname = data[0]
            qlen = int(data[1])
            qstart = int(data[2])
            qqend = int(data[3])
            char = data[4]
            tname = data[5]
            tlen = int(data[6])
            tstart = int(data[7])
            aln_len = int(data[10])
            flag = int(data[11])

            if not flag == 255:

                contig_ref[qname].append(tname)
                contig_ref_aln_length[qname].append(aln_len)
                contig_length[qname] = qlen

        with open(output_folder+output_name, 'w+') as f:
            for k, v in contig_ref.items():
                best = None
                best_len = 0
                c_len = contig_length[k]

                if len(v) > 1:

                    align_sum = 0

                    for ref, l in zip(contig_ref[k], contig_ref_aln_length[k]):

                        align_sum += l

                    if align_sum >= threshold * c_len and best_len < align_sum:
                        best_len = align_sum
                        best = ref

                elif contig_ref_aln_length[k][0] >= threshold * c_len:
                    best = contig_ref[k][0]
                    best_len = contig_ref_aln_length[k][0]
                else:
                    best = 'POOR_MAPPING'
                f.write(f'{k}\t{best}\t{best_len}\n')
