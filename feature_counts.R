#source("http://bioconductor.org/biocLite.R")
#biocLite("Rsubread")
library(Rsubread)

#Feed arguments from shell
args=commandArgs(trailingOnly=TRUE)

#Paired end or not
paired=as.logical(as.numeric(args[2]))
setwd('..')

#pulls list of files to analyze based on extension.
bam_files=list.files(pattern = ("\\.bam$"),recursive=TRUE)

#outputfile name: change "exp_counts.txt" to outputfilename of choice
#output_file="Counts.txt"
output_file="Counts_All_Features.txt"

#use featureCounts to count reads to count transcripts per gene
exp_counts<-featureCounts(files=bam_files,annot.ext=args[1],useMetaFeatures=T,isGTFAnnotationFile=T,GTF.attrType="gene_name",strandSpecific=2,isPairedEnd=paired)

#write counts out in a table
write.table(data.frame(exp_counts$counts,stringsAsFactors=FALSE),file=output_file,quote=FALSE,sep="\t")
