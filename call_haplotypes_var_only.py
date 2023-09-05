#usr/bin/python

#use this to remove variant calls with < 2 alleles, or where the number of variant alleles is < reference alleles
#if  number ref = number variant, assign at random
#in this version if there is only a single allele type of low quality, it is set to null

import glob
import random

#random.seed(100)

file_list = glob.glob("*.vcf")

for fname in file_list:
    datpref = fname[:-4]
    datfile = open(fname, 'r')
    outfile = open(datpref + "_revised", 'w')
    for line in datfile:
        tempvect = line.rstrip('\n').split('\t')
        #use this vector for low quality, i.e. only a single read
        set_to_null = tempvect[0:4]+['.','.','.','.','GT','./.','\n']
        #print(set_to_null)
        #print(tempvect)
        null_line = '\t'.join(set_to_null)
        # if variant called, decide whether to keep it
        # first option: if not nucleotides or high quality, keep GATK's original call
        #print(nuc_counts_vec)
        if tempvect[9] == './.' or tempvect[6] != 'LowQual':
            #print("first test\n")
            newline = line
        elif tempvect[4] == "." and tempvect[6] == 'LowQual':
            newline = null_line
        else:
            #GATK calls variant genotype - reject if based on only a single nucleotide
            temp_7 = tempvect[7].split(';')
            tempv_7_sub = [temp_7[2], temp_7[4]]
            temp_7_new = ';'.join(tempv_7_sub)
            #temp_9 remove extraneous variables from variant call so that e.g. 0:8
            temp_9 = tempvect[9].split(':')
            #print(temp_9)
            two_alleles = temp_9[1].split(',')
            temp_9_sub = ['0', two_alleles[0]]
            temp_9_new = ':'.join(temp_9_sub)
            set_to_ref = tempvect[0:6]+['.', temp_7_new, 'GT:DP:MLPSAC:MLPSAF', temp_9_new, '\n']
            ref_line = '\t'.join(set_to_ref)
            set_to_var = tempvect[0:6]+['.'] + tempvect[7:10] + ['\n']
            var_line = '\t'.join(set_to_var)
            nuc_counts = tempvect[9].split(':')[1]
            nuc_counts_vec = nuc_counts.split(',')
            all_count = tempvect[9].split(':')
            #print(all_count)
            if tempvect[4] == '.':
                ref_line = null_line
            elif (int(nuc_counts_vec[0]) < int(nuc_counts_vec[1]) and int(nuc_counts_vec[1]) > 1):
                newline = '\t'.join(set_to_var)
                # if equal number of variant and reference, generate 0,1 with equal probabilty, randomly assign allele
            elif (int(nuc_counts_vec[0]) == int(nuc_counts_vec[1]) and int(nuc_counts_vec[1]) > 1):
                to_pick = random.randint(0,1)
                #print(to_pick)
                if to_pick == 1:
                    newline = '\t'.join(set_to_var)
                else:
                    newline = '\t'.join(set_to_ref)
                    #print(temp_7)
                    #print(newline)
            elif (int(nuc_counts_vec[0]) > int(nuc_counts_vec[1]) and int(nuc_counts_vec[0]) > 1):
                newline = ref_line
            elif int(all_count[2]) < 2 or (int(nuc_counts_vec[1]) < 2 and int(nuc_counts_vec[0]) < 2):
                newline = null_line
            else:
               newline = line
        outfile.write(newline)
	


