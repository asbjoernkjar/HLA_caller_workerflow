
###### HLA CALLER PIPELINE #######

This pipeline can call HLA types for WES paired end data.
The input is .bam files assembled using hg38. 

The pipeline can call using 3 tools Polysolver, xHLA and Optitype 
It will produce a consensus vote where the HLA-allele with the most votes are chosen. 
If no choice can be made it defaults to using Polysolver


If the alignments are hg19, I recommend using only Polysolver. Since that can run on HG19


###### Setting up the Pipeline #######

The pipeline is controlled by the file config/config.yaml where paths/options should be specified.
Cluster interaction is controlled from the config/slurm/ folder
Cluster.yaml --> needs your cluster profile 
Config.yaml --> controls your snakemake profile (standard options)

Time/memory usage can be specified in the snakemake file or in the command line

### BAMS ###

Feeding BAM files to the program there are two options

Putting the .bam files in results/bam_files
config file should then be fed following path:
	PATH_BAMS_ABSOLUTE: "results/bam_files"

Accessing the files in some other folder, so you don't have to move them

PATH_BAMS_ABSOLUTE: /absolute/path/to/bam_folder

The files will be softlinked into results/bam_files
HOWEVER: singularity needs to be mounted with the origin folder, to acces the files
NOTE: the TRUE absolute path might be different than pwd
	use $readlink -f folder 

### Sample_File ###

The samples are specified by csv file. It should have the following columns:
Sample_ID --> Name of samples, will be used as output names
Bam_file_name --> Name of the bam file for that sample. 
	PATH_BAMS_ABSOLUTE + bam_file_name is used to access files
Race --> race of individual {Caucasian, Black, Asian, Unknown}
	If there is no Race column Unknown is used by default 

NOTE --> the pipeline drops any row with a NA value in any col!!!

###### RUNNING THE PIPELINE ######

Run from conda enviroment with snakemake installed (tested for 6.8.0)
As normal snakemake pipeline:
$snakemake –profile config/slurm -j1

NOTE 
It seems that calling many jobs in parallel on the cluster can result in "communication errors", where temporary files are not created etc. This either result in a crash or a stall.
Reruning the failed jobs (by rerunning snakemake) seems to sovle the problems


###### CHECKING OUTPUT ######

You should check the consensus_call-xxxxxxxxx.out log file. 
This describes if all variants are called via majority vote or if the have defaulted to polysolver
In very rare cases the default polysolver output is not allowed input for netMHCpan 
--> these require manual adjustment to fix. The log file will point to these

###### INSTALLING IMAGES ######

POLYSOLVER
Has been ported to singularity, with a few bug fixes and changes see:
https://github.com/IARCbioinfo/polysolver-singularity
This image can be dowloaded by running 
$singularity pull shub://IARCbioinfo/polysolver-singularity:v4
it should be named polysolver-singularity_v4.sif

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
:-)