# RNAseq
The Gresham Lab RNAseq pipeline

Version 1.0 - 10/18/2017

Usage: The intent of this pipeline is to standardize RNA-seq alignments 
for the Gresham lab and to make the most of the HPC environment.

## To Use:

A. Necessary steps: 
	
	1. Define the SBATCH array size (line 2):
		#SBATCH --array=0-N
		Purpose: SBATCH needs to launch a number of runs equal 
		to the number of fastq files that need to be processed. This 
		command dictates the number of runs to launch. 
		Modify: Change the variable 'N' to the number of fastq files
		contained in the DIR directory. Note the zero indexing, if you have 
		two fastq files, N = 1 and the command is #SBATCH --array=0-1.
	
	2. Define run configuration (line 26-28):
		a. The name of the analysis (line 26):
		analysisName="N"
		Purpose: Name defines output folder name
		b. Is the library Paired End (line 27):
		isPE=N
		Purpose: Paired end and single end libraries require distinct
		operations, this value define the operations applied to your library. 
		Modify: If paired end (PE) set 'N' to 1, if single end (SE) set 'N' to 0.
		c. Does the library have UMIs (line 28):
		isUMI=N
		Purpose: UMI labelled libraries require an additional step to
		to demultiplex, this values dictates whether or not that occurs. 
		Modify: If UMIs were used in library prep, set 'N' to 1, else 
		set 'N' to 0.

	3. Define directory containing fastq files to process (line 34):
		DIR='N'
		Purpose: Directory conatins fastq files to process.
		Modify: Set directory 'N' to the location of your files.
		NOTE: Current sbatch capacities have an upper limit of runs per 
		individual user. If you encouter an 'excessive runs' error you may 
		need to split your files into two seperate directories.

	4. Define compression state (line 56 - 57):
		#FILES=($(ls $DIR/*fastq.gz | perl -pe 's/^.+n0\d_(.+)\.fastq\.gz/$1/g' | sort | uniq))
		FILES=($(ls $DIR/*fastq | perl -pe 's/^.+n0\d_(.+)\.fastq/$1/g' | sort | uniq))
		Purpose: Command collects the filenames of files with a defined 
		file extension, either fastq.gz or fastq.
		Modify: If the files have the gz extension (ie. *.fastq.gz),
		uncomment line 56 and comment line 57. If the files do not have the 
		gz extension (ie. *.fastq), comment line 56 and uncomment line 57.
	
B. Optional Modifications:
	
	1. Reference Genome (line 35):
		hisat_ref='N'
		Purpose: Defines genome to align against.
		Modify: Set 'N' to point to a valid genome. The default genome
		has spike in sequences (BAC, CEL, ERCC) already present.
	
	2. Genome Annotation (line 37):
		gff='N'
		Purpose: Defines annotation file for alignment.
		Modify: Set 'N' to point to a valid 'gff' file. NOTE: Older 
		annotation formats (eg. 'gtf', 'gff2') are not supported. 


##
#	Please submit issues to:
#		https://github.com/GreshamLab/RNAseq/issues
#	--------------------------------------------------- 
##	

## Features for future versions:
	1. General interface python wrapper.
		
