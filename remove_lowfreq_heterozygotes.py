#usr/bin/python

#use this to filter vcfs for spurious "heterozygotes" where the frequency of the second allele is < 25%

import glob
import random

#random.seed(100)

file_list = glob.glob("*.vcf")

for fname in file_list:
    datpref = fname[:-4]
    datfile = open(fname, 'r')
    outfile = open(datpref + "_clean_het", 'w')
    count_het = 0
    count_fix = 0
    for line in datfile:
        tempvect = line.rstrip('\n').split('\t')
        last_element = tempvect[-1]
        if last_element != './.':
            last_element_vec = last_element.split(':')
            gt = last_element_vec[0]
            if gt == '0/1':
                count_het = count_het + 1
                info_element = tempvect[7]
                info_element_vec = info_element.split(';')
                AN = info_element_vec[2]
                DP = info_element_vec[4]
                MQ = info_element_vec[-4]
                MQ0 = info_element_vec[-3]
                revised_pos_7 = ';'.join([AN, DP, MQ, MQ0])
                revised_pos_8 = 'GT:DP'
                nuc_count_vec = last_element_vec[1].split(',')
                n_all_one = int(nuc_count_vec[0])
                n_all_two = int(nuc_count_vec[1])
                total_nucs = n_all_one + n_all_two
                frac_one = n_all_one/(total_nucs + 0.0)
                frac_two = n_all_two/(total_nucs + 0.0)
                #print(nuc_count_vec)
                #print(frac_one)
                #print(frac_two)
                #print("\n")
                if frac_one < 0.25:
                    count_fix = count_fix + 1
                    nuc_count_new = [0, total_nucs]
                    gt_new = '1/1'
                    revised_pos_9 = ':'.join([gt_new, str(total_nucs)])
                    revised_vec = tempvect[0:7] + [revised_pos_7, revised_pos_8, revised_pos_9]
                    new_line = '\t'.join(revised_vec) + '\n'
                elif frac_two < 0.25:
                    nuc_count_new = [total_nucs, 0]
                    gt_new = '0/0'
                    revised_pos_9 = ':'.join([gt_new, str(total_nucs)])
                    revised_vec = tempvect[0:7] + [revised_pos_7, revised_pos_8, revised_pos_9]
                    new_line = '\t'.join(revised_vec) + '\n'
                else:
                    new_line = line
            else:
                new_line = line
        else:
            new_line = line
        outfile.write(new_line)
#print(count_het)
#print(count_fix)
            
        
        
        
    
	


