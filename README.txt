
###### HLA CALLER PIPELINE #######

This pipeline can call HLA types for WES paired end data.
The input is .bam files assembled using hg38. 

The pipeline can call using 3 tools Polysolver, xHLA and Optitype 
It will produce a consensus vote where the HLA-allele with the most votes are chosen. 
If no choice can be made it defaults to using Polysolver


If the alignments are hg19, I recommend using only Polysolver. Since that can run on HG19


###### Setting up the Pipeline #######

The pipeline is controlled by the file config/config.yaml where paths/options should be specified.
Cluster interaction is controlled by your personal snakemake profile
See: https://github.com/Snakemake-Profiles


Time/memory usage can be specified in the snakemake file or in the command line

### BAMS ###

bam paths should not be softlinks 


### Sample_File ###

The samples are specified by csv file. It should have the following columns:
sample_col: 
bam_col: 


RACE_COL --> race of individual {Caucasian, Black, Asian, Unknown}
	If there is no Race column, Unknown is used by default 
    set RACE_COL: null


###### RUNNING THE PIPELINE ######

Run from conda enviroment with snakemake installed (tested for 6.8.0)
As normal snakemake pipeline:
$snakemake –profile your/profile -j1

NOTE 
It seems that calling many jobs in parallel on the cluster can result in "communication errors", where temporary files are not created etc. This either result in a crash or a stall.
Reruning the failed jobs (by rerunning snakemake) seems to sovle the problems


###### CHECKING OUTPUT ######

You should check the consensus_call-xxxxxxxxx.out log file. 
This describes if all variants are called via majority vote or if the have defaulted to polysolver
--> these require manual adjustment to fix. The log file will point to these

###### INSTALLING IMAGES ######

POLYSOLVER
Has been ported to singularity, with a few bug fixes and changes.
using the following build file:

"""
Bootstrap: docker
From: sachet/polysolver:v4
%post
    mkdir /data
    sed -i 's.#!/bin/sh.#!/bin/bash.g' /home/polysolver/scripts/shell_call_hla_type
    sed -i 's.TMP_DIR=\$outDir.TMP_DIR=/tmp.g' /home/polysolver/scripts/shell_call_hla_type
    sed -i 's.TMP_DIR=/home/polysolver.TMP_DIR=/tmp.g' /home/polysolver/scripts/shell_call_hla_type
    sed -i 's.\$SAMTOOLS_DIR./home/polysolver/binaries.g' /home/polysolver/scripts/shell_call_hla_type
    sed -i 's.6:29941260-29945884.chr6:29941260-29945884.g' /home/polysolver/scripts/shell_call_hla_type
    sed -i 's.6:31353872-31357187.chr6:31353872-31357187.g' /home/polysolver/scripts/shell_call_hla_type
    sed -i 's.6:31268749-31272105.chr6:31268749-31272105.g' /home/polysolver/scripts/shell_call_hla_type
    sed -i 's.6:29909037-29913661.chr6:29909037-29913661.g' /home/polysolver/scripts/shell_call_hla_type
    sed -i 's.6:31321649-31324964.chr6:31321649-31324964.g' /home/polysolver/scripts/shell_call_hla_type
    sed -i 's.6:31236526-31239869.chr6:31236526-31239869.g' /home/polysolver/scripts/shell_call_hla_type
"""

run:
$sudo singularity build polysolver-singularity.sif build_file_name
 

xHLA
Has not been ported to singularity, but building it directly from docker works
$singularity build xhla.sif docker://humanlongevity/hla
When making this workflow the image has a core allocation problem, that makes it use
all cores on a node, which results in very large ram usage (>100GB) 
This can be fixed by opening the image and changing line 15 of opt/bin/typing.R: 
options(mc.cores = detectCores()) --> options(mc.cores = 1) 
To open the image run
$singularity build --sandbox xhla/ xhla.sif
Then rebuild using 
$singularity build xhla.sif xhla/

OptiType
Has not been ported to singularity, but building it directly from docker works
$singularity build optitype.sif docker://humanlongevity/hlafred2/optitype


Good Luck :-)
