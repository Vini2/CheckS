<p align="center">
  <img src="CheckS_Logo.png" width="500" title="CheckS logo" alt="CheckS logo">
</p>

# CheckS: Checking the presence of a species in a metagenomics sample with known taxonomic composition

CheckS is a pipeline to check the presence of a species in NGS data (reads and assembled contigs) of a metagenomics sample, provided its reference genome and the reference genomes of the other known species. 

CheckS carries out four analysis steps to determine the presence of a query reference genome as follows.
1. K-mer analysis
2. Coverage analysis
3. Assembly analysis
4. Gene analysis
                                          
## Getting Started

### Dependencies
CheckS requires Python 3.7 (tested on Python 3.7.4). You will need the following tools and python packages installed. Versions tested on are listed as well.
* [Biopython](https://biopython.org/) - version 1.74
* [DSK](https://github.com/GATB/dsk) - version 2.3.3
* [Minimap2](https://github.com/lh3/minimap2) - version 2.18
* [CoverM](https://github.com/wwood/CoverM) - version 0.6.1
* [QUAST](http://bioinf.spbau.ru/quast) - version 5.0.2
* [MetaGeneMark](http://exon.gatech.edu/GeneMark/license_download.cgi)

### Downloading GraphBin2
You can download the latest release of CheckS from [Releases](https://github.com/Vini2/CheckS/releases) or clone the CheckS repository to your machine.

```
git clone https://github.com/Vini2/CheckS.git
```

If you have downloaded a release, you will have to extract the files using the following command.

```
unzip [file_name].zip
```

Now go in to the CheckS folder using the command

```
cd CheckS/
```

### Setting up the environment
We recommend that you use [Conda](https://docs.conda.io/en/latest/) to run CheckS. You can download [Anaconda](https://www.anaconda.com/distribution/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) which contains Conda.

Once you have installed Conda, make sure you are in the `CheckS` folder. Now run the following commands to create a Conda environment and activate it to run CheckS.

```
conda env create -f environment.yml
conda activate checks
```

You will have to download and setup **MetaGeneMark** from [http://exon.gatech.edu/GeneMark/license_download.cgi](http://exon.gatech.edu/GeneMark/license_download.cgi). Next export the path to the MetaGeneMark folder.

```
export MGM_PATH=/path/to/MetaGeneMark_linux_64/mgm
```

Now you are ready to run CheckS.

If you want to switch back to your normal environment, run the following command.

```
conda deactivate
```

## Using CheckS
You can see the usage options of CheckS by typing `./CheckS -h` on the command line. For example,
```
usage: CheckS [-h] --contigs CONTIGS --reads1 READS1 --reads2 READS2
              --query_ref QUERY_REF --other_refs OTHER_REFS --output OUTPUT
              [--ref_ext REF_EXT] [--k K] [--similarity SIMILARITY]
              [--prefix PREFIX] [--nthreads NTHREADS] [-v]

CheckS is a pipeline to check the presence of a species in a metagenomics
sample, provided its reference genome and the reference genomes of the other
known species.

optional arguments:
  -h, --help            show this help message and exit
  --contigs CONTIGS     path to the contigs file
  --reads1 READS1       path to the forward reads file
  --reads2 READS2       path to the reverse reads file
  --query_ref QUERY_REF
                        path to the query reference genome
  --other_refs OTHER_REFS
                        path to the folder with other reference genomes
  --output OUTPUT       path to the output folder
  --ref_ext REF_EXT     extension of the reference genome files. [default:
                        fasta]
  --k K                 k value to run DSK k-mer counting. [default: 25]
  --similarity SIMILARITY
                        similarity threshold for mapping. [default: 0.9]
  --prefix PREFIX       prefix for the output file
  --nthreads NTHREADS   number of threads/cores to use. [default: 8]
  -v, --version         show program's version number and exit
```

## Example Usage

```
./checks --contigs contigs.fasta --reads1 reads_1.fq --reads2 reads_2.fq --query_ref Klebsiella_variicola.fasta --other_refs Reference_Sequences/ --output /path/to/output/
```

## References

[1] DSK: k-mer counting with very low memory usage. Rizk et al. Bioinformatics, Volume 29, Issue 5, 1 March 2013, Pages 652–653.

[2] CoverM - https://github.com/wwood/CoverM

[3] MetaQUAST: evaluation of metagenome assemblies. Mikheenko et al. Bioinformatics, Volume 32, Issue 7, 1 April 2016, Pages 1088–1090.

[4] Minimap2: pairwise alignment for nucleotide sequences. Li. Bioinformatics, Volume 34, Issue 18, 15 September 2018, Pages 3094–3100.

[5] Ab Initio Gene Identification in Metagenomic Sequences. Tang and Borodovsky. In: Nelson K. (eds) Encyclopedia of Metagenomics. Springer, New York, NY.
