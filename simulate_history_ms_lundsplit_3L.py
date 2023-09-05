# generate ms output for analysis of window and site Fst distributions
# vary recombination rate 0, 0.5, 1, 1.5, 2, 3, inf (>3) as input arguments
# new scenario - model population in Lund as a split from Lyon

# include migration in every generation from Lyon to Lund ???
# use Ne estimates for Lund as inferred from Bayesian estimate (Estimate_Pop_Size_Bayesian.R scripts)

import random
import sys
import os
import time
import numpy

start = time.time()

Ne = 827200
theta_site = 0.0077
sample_1800 = 5 #or 5
sample_1933 = 10 #or 14-15
# for poolseq data 3L: 223.595 i.e. 224 raw, sample 142.77 i.e. 143 w/replacement (
# for poolseq data X: 109.66 raw, 70.24 w/replacement (110 and 70)
modern_sample = 240
#nreps = 1000 #for testing
nreps = 100000
tcreate_lund = 0.0009345
tsplit_1933 = 0.0003717 #1230 generations ago, 1230/4Ne
tsplit_1809 = 0.0009339 #3090 generation ago, 3090/4Ne
#make size change slightly later than the split to reduce impact of the split on sampling
t_size_change_1933 = 0.0003710
t_size_change_1809 = 0.0009330
# current estimates 1933 Ne = 2800, 1809 Ne = 1750
# old estimates 1933 Ne = 3200, 1809 Ne = 2400 
size_ratio_1933 = 0.003384913
#size_ratio_1809 = 0.002115571
size_ratio_1809 = 0.003022244
# model population size of 1e+8 or larger to effectively eliminate genetic drift
size_ratio_nodrift = 10
#gene conversion and tract length
gc = 6.25e-8
gc_tract = 518
ploidy = 4 #set to 3 for X chromosome, 4 for autosome

r_rate = sys.argv[1] #recrate = 0, 0.5, 1, 1.5, 2, 3, inf where "inf" > 3
#r_rate = 0.5
recrate_path = "/raid10/mshpak/MuseumFlies/Ref_Files/Rec_Rate/"
output_ms = "demog_ms_output_lundsplit_recrate_"+str(r_rate)
ms_path = '/home/mshpak/To_Install/msdir/'

#placeholders for recombination rate

recrate = []
if ploidy == 4:
   #change _auto_lengths.txt to _Chr3L_lengths.txt (or other chr name)
   # or change to 2R
    with open(recrate_path+'recrate'+str(r_rate)+'_Chr2R_lengths.txt') as recwindows:
        next(recwindows)
        for recraterow in recwindows:
            #print(recraterow)
            recraterow = recraterow.split() # Arm, Start, End, Length, c
            rec_l = recraterow[3:]            
            rec_l[0] = int(rec_l[0]) # len
            rec_l[1] = float(rec_l[1]) # c cM/MB. c = 100 -> 1 Morgan/megabase, 1 crossing-over expected every 1mb
            rec_l.append(rec_l[1]/100000000) # 100000000 -> divide by 100 to convert cM to M. Divide by 1m to convert from Comeron's scale (1MB) to per base.
            rec_l[1] = rec_l[2]*rec_l[0]*Ne*ploidy # multiply by 4N to obtain the population-scaled metric, multiply by rec_l[0] (length) to obtain the locus value
            if rec_l[1] == 0:
                rec_l[2] = ploidy * Ne * gc
            recrate.append(rec_l)
else:
    with open('../../recrate'+r_rate+'_X_lengths.txt') as recwindows:
        next(recwindows)
        for recraterow in recwindows:
            recraterow = recraterow.split() # Arm, Start, End, Length, c
            rec_l = recraterow[3:]
            
            rec_l[0] = int(rec_l[0]) # len
            rec_l[1] = float(rec_l[1]) # c cM/MB. c = 100 -> 1 Morgan/megabase, 1 crossing-over expected every 1mb
            rec_l.append(rec_l[1]/100000000) # 100000000 -> divide by 100 to convert cM to M. Divide by 1m to convert from Comeron's scale (1MB) to per base.
            rec_l[1] = rec_l[2]*rec_l[0]*Ne*ploidy # multiply by 4N to obtain the population-scaled metric, multiply by rec_l[0] (length) to obtain the locus value
            if rec_l[1] == 0:
                rec_l[2] = ploidy * Ne * gc
            recrate.append(rec_l)

#nreps=1            
recrate = random.choices(recrate, k=nreps) # sampling with replacement

#cmd = 'echo "rubbish"'
print(theta_site*rec_l[0])
print('recl[0] = ' + str(rec_l[0]) + '\n')
print('recl[1] = ' + str(rec_l[1]) + '\n')
print('recl[2] = ' + str(rec_l[2]) + '\n')

#replace with demographic multipopulation model with FR as source. Sample FR and its 1809, 1933 splits.
#FR population change to 2400 at 1809, 3250 in 1933
#rec_l[2] (gene conversion rate)
### check if gc is excluded ?????

for rec_i in range(0,nreps):    
    c_size = recrate[rec_i] # randomly obtained from file (with replacement)1
    c_size_ms = str(c_size[1])+' '+str(c_size[0])
    #gene convertion
    if c_size[1] == 0:
        gn_f = c_size[2]
    else:
        gn_f = gc/c_size[2]
    theta =  theta_site * c_size[0] # theta/bp * length    
    gn_tc = gc_tract
    gn_code = str(gn_f)+' '+str(gn_tc)
    #output_ms = 'results_ms'+str(rec_i)+'.txt'
    #use ms 281 1; I -8 0 0 0 0 0 240 14 27
    #original code had -en 0 6 0.360784008 after -I arguments and ' -en ' + str(tsplit_1809) + ' 6 ' + str(size_ratio_1809) + ' -en ' + str(tsplit_1933) + ' 6 ' + str(size_ratio_1933) before end
    cmd = ms_path + 'ms 255 1 -t ' + str(theta) + ' -r ' + str(rec_l[1]) + ' ' + str(rec_l[0]) + str(gc_tract) +  ' -I 9 0 0 0 0 0 0 240 5 10 0 -en 0 1 1.380784346 -en 0 2 0.648736538 -en 0 3 0.554391806 -en 0 4 0.239910374 -en 0 5 0.159309414 -en 0 6 0.360784008 -en 0 7 ' + str(size_ratio_1933) + ' -em 0 1 2 1.77E+01 -em 0 2 1 1.77E+01 -em 0 1 3 1.39E-01 -em 0 3 1 1.39E-01 -em 0 1 4 0 -em 0 4 1 0 -em 0 1 5 0 -em 0 5 1 0 -em 0 1 6 0 -em 0 6 1 0 -em 0 2 3 1.59E+01 -em 0 3 2 1.59E+01 -em 0 2 4 1.26E+01 -em 0 4 2 1.26E+01 -em 0 2 5 1.53E+00 -em 0 5 2 1.53E+00 -em 0 2 6 0 -em 0 6 2 0 -em 0 3 4 2.46E+00 -em 0 4 3 2.46E+00 -em 0 3 5 1.56E+01 -em 0 5 3 1.56E+01 -em 0 3 6 0 -em 0 6 3 0 -em 0 4 5 5.28E+00 -em 0 5 4 5.28E+00 -em 0 4 6 0 -em 0 6 4 0 -em 0 5 6 4.86E+01 -em 0 6 5 4.86E+01 -ej 0.007289253 6 5 -en 0.007289253001 5 0.080531448 -em 0.007289253001 1 2 1.77E+01 -em 0.007289253001 2 1 1.77E+01 -em 0.007289253001 1 3 1.39E-01 -em 0.007289253001 3 1 1.39E-01 -em 0.007289253001 1 4 0 -em 0.007289253001 4 1 0 -em 0.007289253001 1 5 0 -em 0.007289253001 5 1 0 -em 0.007289253001 2 3 1.59E+01 -em 0.007289253001 3 2 1.59E+01 -em 0.007289253001 2 4 1.26E+01 -em 0.007289253001 4 2 1.26E+01 -em 0.007289253001 2 5 1.15E+01 -em 0.007289253001 5 2 1.15E+01 -em 0.007289253001 3 4 2.46E+00 -em 0.007289253001 4 3 2.46E+00 -em 0.007289253001 3 5 2.84E-02 -em 0.007289253001 5 3 2.84E-02 -em 0.007289253001 4 5 2.65E-03 -em 0.007289253001 5 4 2.65E-03 -ej 0.04057006 4 3 -en 0.04057006001 3 0.074500944 -em 0.04057006001 1 2 1.77E+01 -em 0.04057006001 2 1 1.77E+01 -em 0.04057006001 1 3 5.65E-02 -em 0.04057006001 3 1 5.65E-02 -em 0.04057006001 1 5 0 -em 0.04057006001 5 1 0 -em 0.04057006001 2 3 1.10E-01 -em 0.04057006001 3 2 1.10E-01 -em 0.04057006001 2 5 1.15E+01 -em 0.04057006001 5 2 1.15E+01 -em 0.04057006001 3 5 4.87E-02 -em 0.04057006001 5 3 4.87E-02 -ej 0.044395565 5 3 -en 0.044395565001 3 0.595799621 -em 0.044395565001 1 2 1.77E+01 -em 0.044395565001 2 1 1.77E+01 -em 0.044395565001 1 3 5.10E-01 -em 0.044395565001 3 1 5.10E-01 -em 0.044395565001 2 3 1.31E+00 -em 0.044395565001 3 2 1.31E+00 -ej 0.045115152 3 2 -en 0.045115152001 2 0.200985601 -em 0.045115152001 1 2 2.42E+00 -em 0.045115152001 2 1 2.42E+00 -ej 0.046117315 1 2 -en 0.046117315001 2 0.18277622 -en 0.047117966 2 1 -ej ' + str(tcreate_lund) + ' 7 6 -en ' + str(t_size_change_1933) +' 7 ' + str(size_ratio_1809) + ' -ej ' + str(tsplit_1809) + ' 8 7 ' + '-ej ' + str(tsplit_1933) + ' 9 7 ' + '-n 8 ' + str(size_ratio_nodrift) + ' -n 9 ' + str(size_ratio_nodrift) + ' -em 0 7 6 5.00E+01 >>' + output_ms
    
# -en str(t_size_change_1809) 6 0.360784008
# -em ' +str(t_size_change_1933) + ' 9 7 2.40E+00 or 1.20E+00 etc 
# size of Lyon 0.360784008

#cmd = ms_path + 'ms 163 ' + str(nreps) + ' -t ' + str(theta_site*rec_l[0]) + ' -r ' + str(rec_l[1]) + ' ' + str(rec_l[0]) + ' -I 2 27 136 0 -n 2 ' + str(size_ratio) + ' -ej 0.1025 2 1  > ' + output_ms

    os.system(cmd)

end = time.time()
print(end-start)

outfile_name = 'demog_bash_output_lundsplit_' + str(r_rate)
outfile = open(outfile_name, 'w')
outfile.write(cmd+'\n')
