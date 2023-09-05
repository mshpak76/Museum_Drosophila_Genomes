#calculate window fst from ms output

# revised draft takes into account sampling from 3 rather than 2 populations
# also includes twofold poolseq sampling, i.e. effective sample size followed by actual read sampling w/replacement

import re
import sys
import os
import numpy as np
import random
import time
import statistics
import json

start = time.time()

def minor_allele_freq(p):
    if p <= 0.5:
        return(p)
    else:
        return(1-p)
        
        
def major_allele_freq(p):
    if p >= 0.5:
        return(p)
    else:
        return(1-p)


#sum columns of matrix - turn list of lists into numpy object and apply numpy functions        
def colsum(mat):
    arr = np.array(mat)
    csum = np.sum(mat,0)
    return(np.ndarray.tolist(csum))   


# square of every list element
def squarelist(vect):
    return( [i**2 for i in vect])
	

# difference between each element of two lists
def listdiff(v1,v2):
    return( [v1[i] - v2[i] for i in range(len(v1))])


# sum of each element of two lists
def listadd(v1,v2):
    return( [v1[i] + v2[i] for i in range(len(v1))])
    
    
# set to major allele frequencies at all loci
def major_allele_list(v1):
    return( [major_allele_freq(v1[i]) for i in range(len(v1))])
    

#average allele frequency difference
def allele_freq_diff(v1, v2):
    diffv = [abs(v1[i]-v2[i]) for i in range(len(v1))]
    return(statistics.mean(diffv))
    
   

# adjusted reynolds fst for allele (rather than individual) sample
def rey_fst(freq1, freq2, SampleSize1, SampleSize2):
    n1 = SampleSize1
    n2 = SampleSize2
    #sec_allele_1 = 1 - freq1
    #sec_allele_2 = 1 - freq2
    SharedNum = n1*(freq1 - freq1**2) + n2*(freq2 - freq2**2)
    # check if a factor of 1/2 is needed for NumA, as per Reynolds et al 1983
    NumA = (freq1 - freq2)**2
    FracNum = ((n1 + n2)/2)*SharedNum
    FracDen = n1*n2*((n1+n2)/2 - 1)
    frac = FracNum/FracDen
    WholeNum = NumA - frac
    DenFracNum = (n1*n2 - (n1+n2)/2)*SharedNum
    DenFrac = DenFracNum/FracDen
    WholeDen = NumA + DenFrac
    if WholeDen !=0:
        FST = WholeNum/WholeDen
    else:
        FST = 0
    #if (FST > 1):
        #FST = 1
    #elif (FST < 0):
        #FST = 0
    #else:
        #FST = FST
    return([FST, WholeNum, WholeDen])
          
def multilocus_fst(Freqlist1, Freqlist2, SampleSize1, SampleSize2):
    fst_list = [rey_fst(Freqlist1[i], Freqlist2[i], SampleSize1, SampleSize2) for i in range(len(Freqlist1))]
    return(fst_list)
    

#weighted average of multilocus fst    
def weighted_avg_rey_fst(Freqlist1, Freqlist2, SampleSize1, SampleSize2):
    all_dat = multilocus_fst(Freqlist1, Freqlist2, SampleSize1, SampleSize2)
    al_list = [all_dat[i][1] for i in range(len(Freqlist1))]
    albl_list = [all_dat[i][2] for i in range(len(Freqlist1))]
    s_albl = sum(albl_list)
    s_al = sum(al_list)
    if(s_albl > 0):
        temp_ratio = s_al/s_albl
    else:
        temp_ratio = 0
    if temp_ratio < 0:
        temp_ratio = 0
    elif temp_ratio > 1:
        temp_ratio = 1
    else:
        temp_ratio = temp_ratio
    return(temp_ratio)
    
#another approach to multilocus fst - mean of ratios rather than ratio of sums
#def weighted_avg_rey_fst(Freqlist1, Freqlist2, SampleSize1, SampleSize2):
#    all_dat = multilocus_fst(Freqlist1, Freqlist2, SampleSize1, SampleSize2)
#    fst = [all_dat[i][1] for i in range(len(Freqlist1))]
#    #al_list = [all_dat[i][1] for i in range(len(Freqlist1))]
#    #albl_list = [all_dat[i][2] for i in range(len(Freqlist1))]
#    #temp_ratio = sum(al_list)/sum(albl_list)
#    fst_mean = sum(fst)/len(fst)
#    return(fst_mean)
    
    
def singlesite_fst_list(Freqlist1, Freqlist2, SampSize1, SampSize2):
# single locus fst for each site
    fst_temp = [fst(Freqlist1[i], Freqlist2[i], SampSize1, SampSize2)[0] for i in range(len(Freqlist1))]
    numer = [fst(Freqlist1[i], Freqlist2[i], SampSize1, SampSize2)[1] for i in range(len(Freqlist1))]
    denom = [fst(Freqlist1[i], Freqlist2[i], SampSize1, SampSize2)[2] for i in range(len(Freqlist1))]
    return([fst_temp, numer, denom])  

def sample_array_cols(MyMatrix, nelements):
    vmat = []
    TempMat = MyMatrix.T
    #loop over columns rather than rows
    for v in TempMat:
        v = np.ndarray.tolist(v)
        #subv = random.sample(v, nelements)
        subv = random.choices(v, k=nelements)
        vmat = vmat + [subv]
    return(np.array(vmat).T)
    
#def sample_array_cols(matrix, n_result):
#    (n,m) = matrix.shape
#    vmat = numpy.array([n_result, m], dtype= matrix.dtype)
#    for c in range(m):
#        random_indices = numpy.random.randint(0, n, n_result)
 #       vmat[:,c] = matrix[random_indices, c]
 #   return(vmat)
    
    
#####
# analyze ms simulation output, specify parameters    
    
#python ms_fst_reyn.py fst bash_output_0 ms_recrate_output_0
r_rate = sys.argv[1] # 0, 0.5, 1, 1.5, 2, 3, inf

# raw read count = 224, sample size 143
poolseq_sample = 143

input_args = open("demog_bash_output_" + str(r_rate), 'r')
line = next(input_args)
ms_command_list = line.split()
sample_size_1 = int(ms_command_list[15])
sample_size_2 = int(ms_command_list[16])
sample_size_3 = int(ms_command_list[17])
nreps = int(ms_command_list[2])
input_args.close()

#print('sample size 1 = ' + str(sample_size_1) + '\n')
#print('sample size 2 = ' + str(sample_size_2) + '\n')
#print('sample size 3 = ' + str(sample_size_3) + '\n')

# need to revise to compute Fst for populations 1 (modern Lund) vs 3 (1933) and 3 (1933) vs 2 (1809)
# also need to simulate poolseq sampling
#ms_file = open("rubbish", 'r')
ms_file = open("demog_ms_output_recrate_" + str(r_rate), 'r')
start_sim = 0
pos = 0
samp_1 = []
samp_2 = []
samp_3 = []
fst_by_run_2015_1933 = []
fst_by_run_1933_1800 = []
allele_freq_2015_1933 = []
allele_freq_1933_1800 = []
freq_1 = []
freq_3 = []
# samp 1: 2015
# samp 2: 1809
# samp 3: 1933
for line in ms_file:
    #line.strip('\n')
    line = line.rstrip('\n')
    #if line[0:2] == '//':
    if line[0:8] == 'segsites' and int(line.split(':')[1]) > 0:
        start_sim = start_sim + 1
        pos = 0
        if start_sim > 1 and len(samp_2) > 0 and len(samp_1) > 0 and len(samp_3) > 0:
            #print(samp_1)
            #print(samp_2)
            samp_1 = np.array(samp_1)
            samp_1 = sample_array_cols(samp_1, poolseq_sample)
            samp_2 = np.array(samp_2)
            samp_3 = np.array(samp_3)
            vlist_1 = colsum(samp_1)
            vlist_2 = colsum(samp_2)
            vlist_3 = colsum(samp_3)
            vlist_1 = [el1/poolseq_sample for el1 in vlist_1]
            vlist_2 = [el2/sample_size_2 for el2 in vlist_2]
            vlist_3 = [el3/sample_size_3 for el3 in vlist_3]
            meanfreq_1 = sum(vlist_1)/len(vlist_1)
            meanfreq_3 = sum(vlist_3)/len(vlist_3)
            fst_by_run_2015_1933 = fst_by_run_2015_1933 + [weighted_avg_rey_fst(vlist_1, vlist_3, sample_size_1, sample_size_3)]
            fst_by_run_1933_1800 = fst_by_run_1933_1800 + [weighted_avg_rey_fst(vlist_2, vlist_3, sample_size_2, sample_size_3)]
            allele_freq_2015_1933 = allele_freq_2015_1933 + [allele_freq_diff(vlist_1, vlist_3)]
            allele_freq_1933_1800 = allele_freq_1933_1800 + [allele_freq_diff(vlist_2, vlist_3)]
            freq_1 = freq_1 + [meanfreq_1]
            freq_3 = freq_3 + [meanfreq_3]   
        samp_1 = []
        samp_2 = []
        samp_3 = []
    if pos <= (sample_size_1 + 1) and pos > 1 and len(line) > 0 and (line[0] == '0' or line[0]=='1'):
        samp_1 = samp_1 + [[int(x) for x in line]]
    if pos > (sample_size_1 + 1) and pos <= (sample_size_1 +  sample_size_2 + 1) and len(line) > 0 and (line[0] =='0' or line[0]=='1'):
        samp_2 = samp_2 + [[int(x) for x in line]]
    if pos > (sample_size_1 + sample_size_2 + 1) and pos <= (sample_size_1 +  sample_size_2 + sample_size_3 + 1) and len(line)>0 and (line[0]=='0' or line[0]=='1'):
        samp_3 = samp_3 + [[int(x) for x in line]]
    pos = pos + 1
ms_file.close()
    
fst_by_run_2015_1933.sort()
nelements = len(fst_by_run_2015_1933)
mean_fst = sum(fst_by_run_2015_1933)/nelements
mean_af = sum(allele_freq_2015_1933)/nelements
mf1 = sum(freq_1)/nelements
mf3 = sum(freq_3)/nelements

cutoff_2 = [x for x in fst_by_run_2015_1933 if x > 0.2] 
cutoff_3 = [x for x in fst_by_run_2015_1933 if x > 0.3] 
cutoff_4 = [x for x in fst_by_run_2015_1933 if x > 0.4]
cutoff_5 = [x for x in fst_by_run_2015_1933 if x > 0.5]
cutoff_6 = [x for x in fst_by_run_2015_1933 if x > 0.6]
cutoff_7 = [x for x in fst_by_run_2015_1933 if x > 0.7]
cutoff_8 = [x for x in fst_by_run_2015_1933 if x > 0.8]
cutoff_9 = [x for x in fst_by_run_2015_1933 if x > 0.9]

#quantiles for fstcutoffs
p_2 = len(cutoff_2)/nelements
p_3 = len(cutoff_3)/nelements
p_4 = len(cutoff_4)/nelements
p_5 = len(cutoff_5)/nelements
p_6 = len(cutoff_6)/nelements
p_7 = len(cutoff_7)/nelements
p_8 = len(cutoff_8)/nelements
p_9 = len(cutoff_9)/nelements

#afd_13 = allele_freq_diff(vlist_1, vlist_3)
#afd_23 = allele_freq_diff(vlist_2, vlist_3)




#print(fst_by_run_2015_1933)
print("number of elements = " + str(nelements) + '\n')
print("mean fst = " + str(mean_fst) + '\n')
#print("mean allele freq difference " + str(mean_af) + '\n')
#print("mean freq 2015 " + str(mf1) + '\n')
#print("mean freq 1933 " + str(mf3) + '\n')
print("Fst=0.2 " + str(p_2) + " Fst=0.3 " + str(p_3) + " Fst=0.4 " + str(p_4) + " Fst=0.5 " + str(p_5) + " Fst=0.6 " + str(p_6) + " Fst=0.7 " + str(p_7) + " Fst=0.8 " + str(p_8) + " Fst=0.9 " + str(p_9))
#print('\n')

end = time.time()
print("time duration = ")
print(end-start)

outfile=open("fst_distribution_" + str(r_rate), 'w')
outfile.write(str(fst_by_run_2015_1933))

outfile2=open("fst_quantiles_" + str(r_rate), 'w')
outfile2.write(("Fst=0.2 " + str(p_2) + " Fst=0.3 " + str(p_3) + " Fst=0.4 " + str(p_4) + " Fst=0.5 " + str(p_5) + " Fst=0.6 " + str(p_6) + " Fst=0.7 " + str(p_7) + " Fst=0.8 " + str(p_8) + " Fst=0.9 " + str(p_9)))
