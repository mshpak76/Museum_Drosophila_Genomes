#!/usr/bin/python

# # calculate effective sample size for poolseq samples, following Ferretti et al 2013

import sys
import os
import glob

from math import factorial
import sympy

#chr_name = sys.argv[1]
#chr_name = input("Enter name of chromosome: ")
chr_name = "X"

nr_max_40 = 200
nr_max_80 = 400
nr_max_100 = 500

ne_dict_40 = {}
ne_dict_80 = {}
ne_dict_100 = {}

## function to calculate effective sample size as the expectation of unique chromosomal draws
# nr: number reads, nc: number of chromosomes (redo as nc x 2 for diploids = 80)?? Check X chromosome ploidy based on males/females
# possible 24 males/16 females

def expectation_sympy(nr, nc):
    return float(sum(j*factorial(nc)*sympy.functions.combinatorial.numbers.stirling(nr, j)/(factorial(nc - j)*nc**nr) for j in range(1, nc + 1)))
    
def effect_samp_size(nr, nc, nmax):
    if(nr > nmax):
        return expectation_sympy(nmax, nc)
    else:
        return expectation_sympy(nr, nc)
    
for nr40 in range(1,nr_max_40):
    ne_dict_40[nr40] = expectation_sympy(nr40, 40)
    
for nr80 in range(1,nr_max_80):
    ne_dict_80[nr80] = expectation_sympy(nr80, 80)


#for nr100 in range(1, nr_max_100):
#    ne_dict_100[nr100] = expectation_sympy(nr100,100)


#this information isn't used in this script
#poly_sites = []
#poly_file = open("SNP_Sites", 'r')
#next(poly_file)
#for line in poly_file:
#    line = line.rstrip('\n')
#    poly_sites = poly_sites + [line]
    
    
#datfiles = glob.glob('home/mshpak/Lundflies/VCF/Pooled_Round1/CountFile_*_ + chr_name)
datfiles = glob.glob("CountFile_*_" + chr_name)
#datfiles = ['CountFile_SRR8439156_2L']
for fname in datfiles:
    dfile_new = fname + '_new'
    # new file will be based on effective sample size
    dfile = open(fname, 'r')
    newfile = open("Effect_Size_" + fname, 'w')
    #newfile.write("A,C,G,T" + '\n')
    for line in dfile:
        #print("stages1")
        line.rstrip('\n')
        tvec = line.split('\t')
        tvec = list(map(int,tvec))
        rdepth = sum(tvec)
        #effective sample size
        if rdepth==0:
            newfile.write('0' + '\t' + '0' + '\t' + '0' + '\t' + '0'	+ '\n')
        else:
            #print("stages2")
            neff = ne_dict_40[rdepth]
            tfreq = [num/sum(tvec) for num in tvec]
            new_count = [freq*neff for freq in tfreq]
            new_count_convert = [str(x) for x in new_count]
            new_count_str = '\t'.join(new_count_convert)
            newfile.write(new_count_str + '\n')
        
        
