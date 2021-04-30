#!/usr/bin/env python3

import subprocess
import os
import logging

# create logger
logger = logging.getLogger('CheckS 0.1')

def run_dsk(output_folders, k_val, n_cores, output_path):
    for output_folder in output_folders:
        
        logger.debug("Running DSK in "+output_folder)
        
        # Run dsk on ref.fasta
        if not os.path.isfile(output_folder+"/ref.txt"):
            subprocess.run("dsk -nb-cores "+str(n_cores)+" -kmer-size "+str(k_val)+" -file "+output_folder+"/ref.fasta -out-dir "+output_folder, shell=True)
            subprocess.run("dsk2ascii -file "+output_folder+"/ref.h5 -out "+output_folder+"/ref.txt", shell=True)
        
        # Run dsk on refs_combined.fasta
        if not os.path.isfile(output_folder+"/refs_combined.txt"):
            subprocess.run("dsk -nb-cores "+str(n_cores)+" -kmer-size "+str(k_val)+" -file "+output_folder+"/refs_combined.fasta -out-dir "+output_folder, shell=True)
            subprocess.run("dsk2ascii -file "+output_folder+"/refs_combined.h5 -out "+output_folder+"/refs_combined.txt", shell=True)

    if not os.path.isfile(output_path+"/contigs.txt"):
        logger.debug("Running DSK on contigs")
        subprocess.run("dsk -nb-cores "+str(n_cores)+" -kmer-size "+str(k_val)+" -file "+contigs_file+" -out-dir "+output_path, shell=True)
        subprocess.run("dsk2ascii -file "+output_path+"/"+contigs_file_name+".h5 -out "+output_path+"/contigs.txt", shell=True)
    
    if not os.path.isfile(output_path+"/reads.txt"):
        logger.debug("Running DSK on reads")
        subprocess.run("ls -1 "+reads_file_1+" "+reads_file_2+" > "+output_path+"/reads", shell=True)
        subprocess.run("dsk -nb-cores "+str(n_cores)+" -kmer-size "+str(k_val)+" -file "+output_path+"/reads -out-dir "+output_path, shell=True)
        subprocess.run("dsk2ascii -file "+output_path+"/reads.h5 -out "+output_path+"/reads.txt", shell=True)