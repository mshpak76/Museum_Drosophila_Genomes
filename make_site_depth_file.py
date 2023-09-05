#usr/bin/python

# calculate average number of reads per chromosome arm per sample

import sys
import os
import glob
import statistics
import re


# create read depth reference files

#need to revise this to deal with vcf w/o header format, e.g.
  # deal with this format      
 #['3', '2181427', '.', 'A', 'G', '6.02', 'LowQual', 'AC=1;AF=0.500;AN=2;BaseQRankSum=-0.572;DP=8;Dels=0.00;FS=0.000;HaplotypeScore=0.0000;MLEAC=1;MLEAF=0.500;MQ=92.99;MQ0=0;MQRankSum=-1.589;QD=0.75;ReadPosRankSum=-1.589', 'GT:AD:DP:GQ:PL', '0/1:6,2:8:33:33,0,122']
#1       985     .       C       .       .       .       .       GT      ./.
#chr	pos	ref	var

to_rename = ['Yhet', 'mtDNA', '2L', 'X', '3L', '4', '2R', '3R', 'Uextra', '2RHet', '2LHet', '3LHet', '3RHet', 'U', 'XHet', 'Wolbachia']
rename_dict = {'1':'Yhet','2':'mtDNA', '3':'2L', '4':'X', '5':'3L', '6':'4', '7':'2R', '8':'3R','9':'Uextra','10':'2RHet','11':'2LHet','12':'3LHet','13':'3RHet','14':'U','15':'XHet','16':'Wolbachia'}
chr_heading = ' '.join(to_rename)
outfile_prefix = "Persite_Read_Depth_"
#summary_file = open('Site_Read_Depth', 'w')
#os.chdir('home/jeremy/sweden/VCFs/round1')
#for files in glob.glob("*SNPs_filtered.vcf"):
for files in glob.glob("*_shifted.vcf"):
    #print(files)
    files_vect = files.split('.') # remove filename extension
    filename_root = files_vect[0]
    vfile = open(files, 'r')
    #summary_file.write(files+" \n")
    chr_dict = {} # dictionary of chromosomes and depth per site
    #for line in vfile.readlines()[2:]:
    for line in vfile:
        tempvect = line.rstrip('\n').split('\t')
        # pos 0: chromosome, pos 7:info
        if tempvect[4] != '.':
           chrom = tempvect[0]
           info_vect = tempvect[7].split(';')
           depth_part = [x for x in info_vect if re.match('DP=',x)][0]
           #depth = int(info_vect[3].split('=')[1])
           depth = int(depth_part.split('=')[1])
           if chrom in chr_dict:
               chr_dict[chrom]+=[depth]
           else:
               chr_dict[chrom] = [depth]
        else:
           #chr_dict = chr_dict
           chrom = tempvect[0]
           if chrom in chr_dict:
               chr_dict[chrom]+=[0]
           else:
               chr_dict[chrom]=[0]
    chr_summary = {} # dictionary of chromosomes and summary statistics
    #print(chr_summary)
    sample_name = files.split('_')[0] #e.g. H1 etc
    outfile_temp = outfile_prefix + sample_name + "_Chr"
    #for chrom in chr_dict:
    for chrom in ['3','4','5','7','8']:
        outfile_name = outfile_temp + rename_dict[chrom]
        summary_file = open(outfile_name, 'w')
        pos = 1
        for cov_count in chr_dict[chrom]:
            summary_file.write(str(pos) + '\t' + str(cov_count) + '\n')
            pos = pos + 1
            
    	
    	
    	
        #chr_summary[chrom] = [statistics.mean(chr_dict[chrom]), statistics.stdev(chr_dict[chrom]), len(chr_dict[chrom])]
        #summary_file.write("Chromosome "+rename_dict[chrom]+ ": number of Sites = "+str(chr_summary[chrom][2]) +" mean depth = "+str(chr_summary[chrom][0])+" stdev depth = "+str(chr_summary[chrom][1])+"\n") 
    #summary_file.write('\n')
