# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 12:10:31 2017
 sub file. Filters mrna
@author: pieter
"""
import sys

for arg in range(len(sys.argv)):
    if sys.argv[arg] == '-feature':
        feature_name = sys.argv[arg+1]
        feature_name = feature_name.strip()
        
    if sys.argv[arg] == '-gff':
        gff_name = sys.argv[arg+1]
    if sys.argv[arg] == '-i':
        infile_name = sys.argv[arg+1]
    if sys.argv[arg] == '-o':
        outfile_name = sys.argv[arg+1]

    #test line
    #python feature_filter.py -feature CDS -gff /scratch/work/cgsb/genomes/In_house/Fungi/Saccharomyces_cerevisiae/Gresham/GCF_000146045.2_R64_GreshamLabSpikes/GCF_000146045.2_R64_GreshamLabSpikes_genomic.gff -i test_PE_noUMI/Counts.txt -o CDS_feature_counts.txt

gff_file = open(gff_name)
feature_list=[]
all_features=[]

for line in gff_file:
    if line[0]!='#':
        line=line.strip()
        
        check_feature = line.split('\t')[2]
        check_feature = check_feature.strip()
        
        if check_feature not in all_features:
            all_features.append(check_feature)
        
        if check_feature==feature_name:
            deets = line.split('\t')[8]

            if 'gene=' in deets:
                check_feature=deets.split('gene=')[1].split(';')[0]
            else:
                check_feature=deets.split('ID=')[1]
                check_feature=deets.split(';')[0]
            if len(check_feature) > 0:
                feature_list.append(check_feature)
                               
gff_file.close()

if len(feature_list)<10:
    print('Only '+str(len(feature_list))+' features found. Try one of these features, instead:\n')
    print(all_features)

infile = open(infile_name)
outfile = open(outfile_name,'w')

for line in infile:
    if (line[0:4] == 'BAC-') or (line[0:4] == 'CEL-') or (line[0:5] == 'ERCC-') or ('.bam' in line):
        outfile.write(line)
    else:
        check_feature = line.split('\t')[0]
        if check_feature in feature_list:
            outfile.write(line)

infile.close()
outfile.close()
