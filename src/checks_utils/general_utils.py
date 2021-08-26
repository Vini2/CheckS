#!/usr/bin/env python3

import subprocess
import os


def get_ref_names(ref_files, ext):

    ref_names = []

    for file in ref_files:
        ref_names.append(file.split("/")[-1][:-(len(ext)+1)])

    return ref_names


def make_output_folders(output_path, ref_names):

    output_folders = []

    for ref in ref_names:
        output_folders.append(output_path+ref+"/")
        subprocess.run("mkdir -p "+output_path+ref, shell=True)

    return output_folders


def prep_files(ref_files, output_folders):

    for i in range(len(ref_files)):
        refs_excluding_reference = [x for x in ref_files if x != ref_files[i]]
        subprocess.run("cat "+" ".join(refs_excluding_reference) +
                       " > "+output_folders[i]+"/refs_combined.fasta", shell=True)
        subprocess.run("cp "+ref_files[i]+" " +
                       output_folders[i]+"/ref.fasta", shell=True)
