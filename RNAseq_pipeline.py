# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 16:36:19 2017

This script acts as a wrapper to generate a bash file with certain parameters 
based on user input. By defualt SBATCH is automatically launched using the 
generated command file.

@author: Pieter Spealman. email pspealman@nyu.edu
"""

### import section. sys handles user IO, subprocess makes system calls.
import sys
import subprocess

### defaults section. Sets default parameters for wrapper.
input_dir = ''
output_dir = ''
isPE = 0
isUMI = 0
file_number = 1
file_ext='fastq'

genome_name='/scratch/work/cgsb/genomes/In_house/Fungi/Saccharomyces_cerevisiae/Gresham/GCF_000146045.2_R64_GreshamLabSpikes/GCF_000146045.2_R64_GreshamLabSpikes_genomic_hisat_index/genome'
gff_name='/scratch/work/cgsb/genomes/In_house/Fungi/Saccharomyces_cerevisiae/Gresham/GCF_000146045.2_R64_GreshamLabSpikes/GCF_000146045.2_R64_GreshamLabSpikes_genomic.gff'
polyG = False
blank_lines = False
no_launch = False

#function handles help commands
def help_dialog():
    print('\n')
    print('#=============================================================#')
    print('\tThis script takes input to generate a standard shell script')
    print('that aligns RNAseq data and performs QC. The default is assumed')
    print('to be single end (unpaired), fastq files, without UMIs.\n')
    print('#=============================================================#')
    print('Usage:')
    print('python RNAseq_pipeline.py [-pe] [-umi] [-gzip] -n file_number -i input_directory -o output_directory\n')
    print('\t\tCommands:')
    print('-n : Number of files (.fastq or fastq.gzip), present in directory.')
    print('\t this number is used to launch arrays. Begin counting at one.')
    print('-i : Input directory which holds the files to be processed')
    print('-o : Output directory will hold the aligned files and logs.\n')
    print('\t\tOptional Commands:')
    print('-pe : Paired end read library. Default is unpaired reads.')
    print('-umi : UMI adapter library. UMItools will be run on processed reads.')
    print('\t\tTurned off by defualt.')
    print('-gzip : FASTQ files have been compressed with gzip. By defualt files')
    print('\t\tare asumed to be uncompressed. Other compression formats not supported\n')
    print('\t\tExtended optional commands.')
    print('[-polyG]')
    print("-polyG : Nextera's two color chemistry reads the absence of signal as 'G'")
    print("\t\tas such strings of 'GGG GGG' are artificially inflated. This option removes")
    print("\t\tthose cases which have 'GGG GGG' UMI labels. Usage implies -umi.")
    print('[-blank_lines]')
    print('-blank_lines : Older seqeuncing files conatin an erroneous new line')
    print('\t\tcharacter, while this bug has been fixed, this option allows for')
    print('\t\tbackwards compatibility.')
    print('[-gff gff_name]')
    print("-gff gff_name: 'gff_name' specifies the custom gff reference file to use.")
    print('[-genome genome_name]')
    print("-genome genome_name: genome_name specifies the custom reference genome to use.")
    print('[-no_launch]')
    print('-no_launch : This option does not invoke SBATCH after generating the')
    print('\t\tcommand file.')
    print('')    
    print('test se umi:\n\t python RNAseq_pipeline.py -umi -gzip -n 2 -i /scratch/cgsb/gencore/out/Gresham/2015-04-24_C5726ACXX/umi/8 -o test_SE_UMI')
    print('test pe no umi:\n\t python RNAseq_pipeline.py -pe -gzip -n 2 -i /scratch/cgsb/gresham/LABSHARE/Scripts/Pipelines/testdata_PE -o test_PE_noUMI')
    sys.exit('')    

### user input section
# parse sys.argv for help requests
if ('-h' in sys.argv) or ('-help' in sys.argv) or (len(sys.argv)<1) or ('-man' in sys.argv) or ('?' in sys.argv):
    help_dialog()

# parse sys.argv for commands and options
for arg in range(len(sys.argv)):
    if sys.argv[arg] == '-n': 
        file_number = int(sys.argv[arg+1])
    if sys.argv[arg] == '-i':
        input_dir = sys.argv[arg+1]
    if sys.argv[arg] == '-o':
        output_dir = sys.argv[arg+1]
        
    if (sys.argv[arg] == '-pe') or (sys.argv[arg] == '--paired_end') :
        isPE = 1
    if (sys.argv[arg] == '-umi') or (sys.argv[arg] == '--umi') :
        isUMI = 1
    if (sys.argv[arg] == '-gzip'):
        file_ext='fastq.gz'
        
    if sys.argv[arg] == '-gff':
        gff_name = sys.argv[arg+1]
    if sys.argv[arg] == '-genome':
        genome_name = sys.argv[arg+1]
    if sys.argv[arg] == '-no_launch':
        no_launch = True
    if sys.argv[arg] == '-polyG':
        isUMI = 1
        polyG = True
    if sys.argv[arg] == '-blank_lines':
        blank_lines = True

# eval input for errors
if file_number > 1:
    array_start = 0
    array_stop = file_number-1       

if input_dir =='':
    print("Input directory error. Please specify a input directory using the '-i' argument")
    help_dialog()
    
if output_dir =='':
    print("Output directory error. Please specify an output directory using the '-o' argument")
    help_dialog()
    
outfile=open('rnaseq-array-pipeline.sh','w')
analysisName = output_dir

# function generates module load section with difference based on UMI 
def module_load():
    outline=('module purge\n'+
    'module load cutadapt/intel/1.12\n'+
    'module load hisat2/intel/2.0.5\n'+
    'module load samtools/intel/1.3.1\n'+
    'module load r/intel/3.3.2\n'+
    'module load qualimap/2.2.1\n'+
    'module load gffread/intel/0.9.8c\n'+
    'module load biopython/intel/1.70\n')
    outfile.write(outline)
    
    if isUMI==0:
        outline =('module load trim_galore/0.4.4\n\n')
        outfile.write(outline)

    if isUMI==1:
        outline=('module load fastqc/0.11.5\n'+
        'module load umi_tools/0.4.2\n\n')
        outfile.write(outline)

# function generates setup section with difference based on array_size (n) and file_ext (-gzip)
# This function is always called.
def write_setup_section():
    
    outline='#!/bin/bash\n'
    outfile.write(outline)

    if array_stop == 0:
        outline = ('#SBATCH --array=%s ##specify number of samples\n')%(array_start)
        outfile.write(outline)   

    if array_stop != 0:
        outline = ('#SBATCH --array=%s-%s ##specify number of samples\n')%(array_start,array_stop)
        outfile.write(outline)
    
    outline = ('#SBATCH --output=%A_%a.o\n'+
    '#SBATCH --error=%A_%a.e\n'+
    '#SBATCH --job-name=%a_rnaseq-array-pipeline\n'+
    '#SBATCH --nodes=1\n'+
    '#SBATCH --ntasks-per-node=1\n'+
    '#SBATCH --mem=60G\n'+
    '#SBATCH --time=04:00:00\n'+
    '#SBATCH --mail-user=${USER}@nyu.edu\n'+
    '#SBATCH --mail-type=BEGIN,END,FAIL\n\n')
    outfile.write(outline)
    
    module_load()

    outline = ('#Configuration\n'+
    'analysisName=%s\n'+
    'isPE=%s #0 if SE\n'+
    'isUMI=%s #1 if UMIs used\n\n')%(analysisName,isPE,isUMI)
    outfile.write(outline)
    
    DIR=input_dir
    hisat_ref=genome_name
    gff=gff_name
    
    outline=('DIR=%s\n'+
    'hisat_ref=%s\n'+
    'gff=%s\n\n')%(DIR,hisat_ref,gff)
    outfile.write(outline)
    
    outline=('#Make a local directory and out-of-stream resource copies of this script for records\n'+
    'mkdir -p ${analysisName}\n'+
    'mkdir -p ${analysisName}/logs\n'+
    'cp rnaseq-array-pipeline.sh "${analysisName}/logs/rnaseq-array-pipeline.$(date +%F_%R).sh"\n'+
    'cd ${analysisName}\n'+
    '\necho "Generating the gtf file from the gff, using gffread"\n'+
    'time gffread -E ${gff} -T -o- > temp.gtf\n'+
    'gtf1="temp.gtf"\n\n')
    outfile.write(outline)
    
    if file_ext=='fastq.gz':
        outline=("FILES=($(ls $DIR/*fastq.gz | perl -pe 's/^.+n0\d_(.+)\.fastq\.gz/$1/g' | sort | uniq))\n")
    
    if file_ext=='fastq':
        outline=("FILES=($(ls $DIR/*fastq | perl -pe 's/^.+n0\d_(.+)\.fastq/$1/g' | sort | uniq))\n\n")
    
    outfile.write(outline)
    
    outline=('## $ID will hold the sample_id\n'+
    'ID=${FILES[$SLURM_ARRAY_TASK_ID]}\n\n')
    outfile.write(outline)
    
    outline=('echo "This is the FILE Im running"\n'+
    'echo ${FILES}\n'+
    'echo "This is its ID"\n'+
    'echo ${ID}\n'+
    'echo "This is the ARRAY number Ive been assigned"\n'+
    'echo ${SLURM_ARRAY_TASK_ID}\n'+
    'echo "Im making this directory and jumping into it: ${ID}"\n'+
    'mkdir $ID\n'+
    'cd $ID\n'+
    '\n#Update gtf variable as it only exists once and is in analysis directory\n'+
    'gtf="../"${gtf1}\n'+
    'echo "I think my GTF is at ${gtf}"\n\n')    
    outfile.write(outline)

# function generates umi_tools section, only called with -umi option
def umi_tools_section():
    
    outline = ('\necho "You want me to dedup? Okay, deduping with umi_tools, then sorting and deleting the unsorted bam."\n'+
    'time umi_tools dedup --umi-separator=: --method=directional --output-stats=${ID}.stats -I ${ID}.bam -S ${ID}.dedup.bam -L ${ID}.dedup.log\n'+
    'time samtools sort -o ${ID}.sorted.dedup.bam ${ID}.dedup.bam\n'+
    'rm ${ID}.dedup.bam\n\n')
    outfile.write(outline)
    
    outline = ('##Determine counts\n'+                                                                                                                                      
    'echo "There were "$( samtools view ${ID}.bam | wc -l )" counts before umi-deduplication"\n'+
    'echo "There are "$( samtools view ${ID}.sorted.dedup.bam | wc -l )" counts after umi-deduplication"\n\n')
    outfile.write(outline)
    
    outline=('\necho "Deleting a lot of intermediate files, indexing."\n'+
    'rm ${ID}.bam\n'+
    'rm ${ID}.bam.bai\n'+
    'mv ${ID}.sorted.dedup.bam ${ID}.bam\n'+
    'time samtools index ${ID}.bam\n\n')
    outfile.write(outline)

# function generates body, only called with -pe and no -umi
# pe invokes trim_galore as no umi adapters are needed, 
# trim_galore autoidentifies Nextera adapters for trimming. 
def paired_end_no_umi():
    
    outline=('## full paths to paired mate files\n'+
    'INPUT_1=$(ls $DIR/*n01_$ID.%s)\n'+
    'INPUT_2=$(ls $DIR/*n02_$ID.%s)\n\n')%(file_ext,file_ext)
    outfile.write(outline)

    if blank_lines:
            outline = ("awk 'NF' $INPUT_1 > INPUT_1.fastq\n"+
            "INPUT_1=$(ls INPUT_1.fastq)\n"+
            "awk 'NF' $INPUT_2 > INPUT_2.fastq\n"+
            "INPUT_2=$(ls INPUT_2.fastq)\n\n")
            outfile.write(outline)
    
    outline=('\necho "Now Im trimming the reads for adapters, quality?"\n'+
    'time trim_galore --paired --fastqc $INPUT_1 $INPUT_2\n'+
    'time hisat2 -p ${SLURM_CPUS_ON_NODE} -q -x ${hisat_ref} -1 *n01_${ID}_val_1.fq.gz -2 *n02_${ID}_val_2.fq.gz -S ${ID}.sam\n\n')
    outfile.write(outline)
        
    outline=('\necho "Okay, now I sort and index that bam"\n'+
    'time samtools sort -o ${ID}.bam ${ID}.sam\n'+
    'time samtools index ${ID}.bam\n'+
    '\necho "And delete that file, and fastqz files as well."\n'+
    'rm ${ID}.sam\n'+
    'rm *n01_${ID}_val_1.fq.gz\nrm *n02_${ID}_val_2.fq.gz\n\n')       
    outfile.write(outline)

    outline=('\necho "Now doing qualimaps"\n'+
    'unset DISPLAY\n'+
    'time qualimap bamqc -bam ${ID}.bam -gff ${gff} -c -outdir ${ID}_stats -outfile ${ID}\n'+
    'time qualimap rnaseq -bam ${ID}.bam -gtf ${gtf} -pe -outdir rnaseq_qc_results -outfile ${ID}\n\n')
    outfile.write(outline)

# function generates body, only called with -pe and -umi
# because we are using umi we need to identify adapters, for this we call 'scrape_adapter'
# and pass the results to cutadapt. After alignment we call umi_tools.
def paired_end_with_umi():
    
    outline=('## full paths to paired mate files\n'+
    'INPUT_1=$(ls $DIR/*n01_$ID.%s)\n'+
    'INPUT_2=$(ls $DIR/*n02_$ID.%s)\n\n')%(file_ext,file_ext)
    outfile.write(outline)
    
    if blank_lines:
            outline = ("awk 'NF' $INPUT_1 > INPUT_1.fastq\n"+
            "INPUT_1=$(ls INPUT_1.fastq)\n"+
            "awk 'NF' $INPUT_2 > INPUT_2.fastq\n"+
            "INPUT_2=$(ls INPUT_2.fastq)\n\n")
            outfile.write(outline)
            
    if polyG:
        outline = ('time python ../../scrape_adapter.py -polyG -i $INPUT_1 -o Read1_filter_polyG.fastq\n'+
        'time python ../../scrape_adapter.py -polyG -i $INPUT_2 -o Read2_filter_polyG.fastq\n'+
       'echo "I assume youre using the trumiseq adapters, so trimming"\n'+
        'time python ../../scrape_adapter.py -u -i Read1_filter_polyG.fastq -o Read1_Adapter.fa\n'+
        'time python ../../scrape_adapter.py -p -i Read2_filter_polyG.fastq -o Read2_Adapter.fa\n'+
        'time cutadapt -a file:Read1_Adapter.fa -A file:Read2_Adapter.fa -o ${ID}_1.trim -p ${ID}_2.trim $INPUT_1 $INPUT_2\n'+
        'time fastqc -o ${ID}_postTrimming -f fastq ${ID}_1.trim ${ID}_2.trim\n'+
        'time hisat2 -p ${SLURM_CPUS_ON_NODE} -q -x ${hisat_ref} -1 ${ID}_1.trim -2 ${ID}_2.trim -S $ID.sam\n\n')
        outfile.write(outline)
            
    if not polyG:
        outline = ('echo "I assume youre using the trumiseq adapters, so trimming"\n'+
        'time python ../../scrape_adapter.py -u -i $INPUT_1 -o Read1_Adapter.fa\n'+
        'time python ../../scrape_adapter.py -p -i $INPUT_2 -o Read2_Adapter.fa\n'+
        'time cutadapt -a file:Read1_Adapter.fa -A file:Read2_Adapter.fa -o ${ID}_1.trim -p ${ID}_2.trim $INPUT_1 $INPUT_2\n'+
        'time fastqc -o ${ID}_postTrimming -f fastq ${ID}_1.trim ${ID}_2.trim\n'+
        'time hisat2 -p ${SLURM_CPUS_ON_NODE} -q -x ${hisat_ref} -1 ${ID}_1.trim -2 ${ID}_2.trim -S $ID.sam\n\n')
        outfile.write(outline)

    outline = ('\necho "Okay, now I sort and index that bam"\n'+
    'time samtools sort -o ${ID}.bam ${ID}.sam\n'+
    'time samtools index ${ID}.bam)\n'+
    '\necho "And delete that file, and fastqz files as well."\n'+
    'rm ${ID}.sam\n'+
    'rm ${ID}_1.trim\n'+
    'rm ${ID}_2.trim\n\n')
    outfile.write(outline)
    
    umi_tools_section()

# function generates body, only called with no -pe and no -umi (default)
# invokes trim_galore as no umi adapters are needed, 
# trim_galore autoidentifies Nextera adapters for trimming.   
def single_end_no_umi():
    
    outline = ('INPUT_1=$(ls $DIR/*n01_$ID.%s)\n\n')%(file_ext)
    outfile.write(outline)
    
    if blank_lines:
        outline = ("awk 'NF' $INPUT_1 > INPUT_1.fastq\n"+
        "INPUT_1=$(ls INPUT_1.fastq)\n\n")
        outfile.write(outline)
    
    outline = ('\necho "Trimgaloring that ${INPUT_1}, then aligning"\n'+
    'time trim_galore --fastqc $INPUT_1\n'+
    'time hisat2 -q -x ${hisat_ref} -1 $INPUT_1 $ID.sam\n'+
    '\necho "sort and index"\n'+
    'time samtools sort -o ${ID}.bam ${ID}.sam\n'+
    'time samtools index ${ID}.bam\n'+
    '#Clean up .sam and .fastq files\n'+
    'rm ${ID}.sam\n\n')
    outfile.write(outline)

# function generates body, only called with no -pe and -umi
# because we are using umi we need to identify adapters, for this we call 'scrape_adapter'
# and pass the results to cutadapt. After alignment we call umi_tools.
def single_end_with_umi():
    
    outline = ('INPUT_1=$(ls $DIR/*n01_$ID.%s)\n\n')%(file_ext)
    outfile.write(outline)

    if blank_lines:
        outline = ("awk 'NF' $INPUT_1 > INPUT_1.fastq\n"+
        "INPUT_1=$(ls INPUT_1.fastq)\n\n")
        outfile.write(outline)    
    
    if polyG:
        outline = ('time python ../../scrape_adapter.py -polyG -i $INPUT_1 -o Read1_filter_polyG.fastq\n'+
        '\necho "Trimming that ${INPUT_1} for adapters, quality? aligning"\n'+
        '\necho ${ID}_1.trim\n'+
        'time python ../../scrape_adapter.py -u -i Read1_filter_polyG.fastq -o Read1_Adapter.fa\n'+
        'time cutadapt -a file:Read1_Adapter.fa -o ${ID}_1.trim $INPUT_1\n'+
        'time hisat2 -q -x ${hisat_ref} -U ${ID}_1.trim -S $ID.sam\n\n')
        outfile.write(outline)
        
    if not polyG:
        outline = ('\necho "Trimming that ${INPUT_1} for adapters, quality? aligning"\n'+
        '\necho ${ID}_1.trim\n'+
        'time python ../../scrape_adapter.py -u -i $INPUT_1 -o Read1_Adapter.fa\n'+
        'time cutadapt -a file:Read1_Adapter.fa -o ${ID}_1.trim $INPUT_1\n'+
        'time hisat2 -q -x ${hisat_ref} -U ${ID}_1.trim -S $ID.sam\n\n')
        outfile.write(outline)

    outline = ('echo "sort and index"\n'+
    'time samtools sort -o ${ID}.bam ${ID}.sam\n'+
    'time samtools index ${ID}.bam\n'+
    '#Clean up .sam and .fastq files\n\n'+
    'rm ${ID}.sam\n')
    outfile.write(outline)
    
    umi_tools_section()    

# function generates rsubread section. This is always called.
def rsubread():
    outline=('############################\n# Generate table of counts #\n############################\n#Call Rsubread\n'+
    '\necho "and doing rsubread"\n'+
    'time Rscript --vanilla ../../feature_counts.R ${gtf1} ${isPE} > ../logs/rsubread.log\n'+
    'time python ../../feature_filter.py -feature CDS -gff ${gff} -i ../Counts_All_Features.txt -o ../Counts_CDS_Features.txt\n'+
    '\necho "Completed now moving the STDOUT and STDERR logs from sbatch into the above directory."\n\n')
    outfile.write(outline)

# function generates clean-up section. This is always called.
def clean_up():
    outline=('#################################################\n# Move output and error files to logs directory #\n'+
    '#################################################\n'+
    'echo $(pwd)\echon $(ls)\n'+
    'mv ../../${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.o ../logs\n'+
    'mv ../../${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.e ../logs\n'+
    '\nexit 0;')
    outfile.write(outline)
    
### Function call section. Takes the parameters together to decide output generated.
write_setup_section()

if isPE==1: 
    if isUMI==0:
        paired_end_no_umi()
    if isUMI==1:
        paired_end_with_umi()

if isPE==0:
    if isUMI==0:
        single_end_no_umi()
    if isUMI==1:
        single_end_with_umi()
    
rsubread()
clean_up()
outfile.close()

if not no_launch:
    output = subprocess.check_output(['sbatch rnaseq-array-pipeline.sh'],stderr=subprocess.STDOUT,shell=True)
    print(output)

if no_launch:
    print('-no_launch option suppressed SBATCH launch. Command file generated.')
