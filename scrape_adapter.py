# -*- coding: utf-8 -*-
"""
Created on Tue Oct 17 13:11:16 2017

The purpose of this script is to scrape the fastq file and enter the sequencer 
defined adapter to a list. We then use this adapter list to generate a '.fa' 
file that cutadapt can use to perform trimming.

@author: pieter
"""
import sys
import subprocess
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

#set defaults
paired_arg = '-u'
adapter_list=[]
polyG = False

#collect arguments
if (sys.argv[1] == '-h') or (sys.argv[1] == '-help'):
    print('scrape_adapter is designed to be used with fastq files.\nThese can be either in uncompressed of gzipped (".gz") formats.')
    print('Invoke as:\n\tpython scrape_adapter.py [-polyG] [-p] -i fastq_file_name -o output_file_name')
    print('-------------------------------------------------------------------')
    print('Optional Commands:')
    print('-p : paired end reads')
    print('-polyG : Alternative Runmode. Removes all reads with a polyG UMI.')
    print('\t\tNo adapters are counted when using -polyG run mode.')
    print('')
    exit()

for arg in range(len(sys.argv)-1):
    if sys.argv[arg] == '-p':
        paired_arg = '-p'  
    if sys.argv[arg] == '-i':
        infile_name = sys.argv[int(arg+1)]
    if sys.argv[arg] == '-o':
        outfile_name = sys.argv[int(arg+1)]
    if sys.argv[arg] == '-polyG':
        polyG = True

#test to see if input is gzipped
#if it is use zcat to make temp uncompressed file
if infile_name.rsplit('.',1)[1]=='gz':
    is_gzipped=True
else:
    is_gzipped=False

# process gzipped file then open temp
if is_gzipped:
    bashCommand = ('zcat %s > temp.fastq')%(infile_name)
    output = subprocess.check_output([bashCommand],stderr=subprocess.STDOUT,shell=True)   
    infile_o = ('temp.fastq')

#open uncompressed file
if not is_gzipped:
    infile_o = (infile_name)

# search for and remove instances where UMI is 'GGGGGG'
def remove_polyG(infile_o):
    infile = open(infile_o)
    
    
    filter_polyG = 0
    polyG_ct = 0
    
    for line in infile:

        if filter_polyG != 0:
            if filter_polyG == 3:
                filter_polyG = 0
                print(line)
                print(filter_polyG)
            if (filter_polyG > 0) and (filter_polyG < 3) :
                filter_polyG+=1
                print(line)
                print(filter_polyG)
                
        if filter_polyG == 0:
            if (line[0]=='@'):
                s_line = line.strip()
                adapter = s_line.rsplit(':',1)[1]
                if ('GGGGGG' in adapter):
                    filter_polyG = 1
                    polyG_ct+=1
                    print(line)
                    print(filter_polyG)
                    
                if ('GGGGGG' not in adapter):
                    outfile.write(line)
                    
            if (line[0]!='@') and (filter_polyG == 0):
                    outfile.write(line)
                

    infile.close()
    print('PolyG count: '+str(polyG_ct))

# parse file keeping only highest frequency adapter
def scrape_file(infile_o):
    infile = open(infile_o)
    
    adapter_dict = {}
    high_score=0
    winner=''
    #parse file
    for line in infile:
        if line[0]=='@':
            line = line.strip()
            adapter = line.rsplit(':',1)[1]
            if ('N' not in adapter):
                if (adapter in adapter_dict): 
                    adapter_dict[adapter]+=1
                if (adapter not in adapter_dict):
                    adapter_dict[adapter]=1
    print('Scraped the following adapters these many times:')
    print(adapter_dict)
    
    for key, value in adapter_dict.items():
        if value > high_score:
            high_score=value
            winner=key
            
    adapter_list.append(winner)
    infile.close()
    return(adapter_list)

#write out adapter        
def write_outadapter(adapter_list,runmode):
    if runmode =='fwd':
        count = 0
        print('Using Adapter' +str(adapter_list))
        for adapter in adapter_list:
            revadapter = Seq(adapter, IUPAC.unambiguous_dna)
            adapter = revadapter.reverse_complement()
            outline=('>FwdAdapter_%s\na%sagatcggaagagcacacgtctgaactccagtcacNNNNNNatctcgtatgccgtcttctgcttg\n')%(count,adapter)
            print('Trimming with:')
            print(outline)
            count+=1
            outfile.write(outline)
            
    if runmode =='rev':
        count = 0
        print('Using Adapter' +str(adapter_list))
        for adapter in adapter_list:
            adapter = Seq(adapter, IUPAC.unambiguous_dna)
            #adapter = revadapter.reverse_complement()
            outline=('>RevAdapter_%s\naatgatacggcgaccaccgagatctacactctttccctacacgacgctcttccgatct%sa\n')%(count,adapter)
            print('Trimming with:')
            print(outline)
            count+=1
            outfile.write(outline)

#actual command statements
if not polyG:
    if paired_arg=='-u':
        adapter_list = scrape_file(infile_o)
        
        if adapter_list:
            outfile = open(outfile_name,'w')
            write_outadapter(adapter_list,'fwd')
    
    if paired_arg=='-p':
        adapter_list = scrape_file(infile_o)
            
        if adapter_list:
            outfile = open(outfile_name,'w')
            write_outadapter(adapter_list,'rev')

if polyG:
   outfile = open(outfile_name,'w')
   remove_polyG(infile_o)

#clean up
outfile.close()
if is_gzipped:
    output = subprocess.check_output(['rm temp.fastq'],stderr=subprocess.STDOUT,shell=True)
