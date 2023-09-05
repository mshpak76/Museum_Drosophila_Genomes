# calculate difference in allele frequency expected based on model effective population size
# at 1933-2015 time intervals

# allow singletons (using high quality sequences with masked IBD regions only)
# redo for Chr 2L, 2R, 3L, 3R exclude low recombination regions. Priroitize 3L and 2R, i.e. fewer inversions

#mutation rate per nucleotide (Huang et al 2020)
dm_mut = 5.21*10**(-9)
#Ne = 1000

#equilibrium allele frequency
Mutation_Drift = function(freq, mu, ne){
	pop_mut = 2*ne*mu
	density_fun = (freq*(1-freq))**(pop_mut-1)
	return(density_fun)
}


Freq_Interval = c(1:99)/100

mut_drift_test = function(f1){
	return(Mutation_Drift(f1, dm_mut, Ne))
	}

#list arguments that are not part of apply after function
#Ne=10000
#Mut_Drift_Dist = sapply(Freq_Interval, Mutation_Drift, dm_mut, Ne)
#Mut_Drift_Dist = Mut_Drift_Dist/sum(Mut_Drift_Dist)

#Uni_Dist = rep(1,99)/99

#Mut_Drift_Dist = c(0, Mut_Drift_Dist, 0)

#Bayesian estimate of actual allele frequency
#freq_pop: population frequency, prior_freq_pop is prior probability
Bayes_Allele_Freq = function(posit, prior_freq_pop, freq_observed, sample_size){
	nobs = round(freq_observed*sample_size)
	nlen = length(prior_freq_pop)
	condprob = dbinom(nobs, prob = posit/(nlen+1), size = sample_size)
	return(condprob*prior_freq_pop[posit])
}

#fobs = .6
#n_alleles = 26
#bayes_intervals = c(1:99)
#Estimate_Allele_Freq_Dist = sapply(bayes_intervals, Bayes_Allele_Freq, Mut_Drift_Dist, fobs, n_alleles)
#Estimate_Allele_Freq_Dist = Estimate_Allele_Freq_Dist/sum(Estimate_Allele_Freq_Dist)
#Cumulative_Est_Allele_Freq = cumsum(Estimate_Allele_Freq_Dist)

#length(which(Cumulative_Est_Allele_Freq <= runif(1,0,1))/100
#Allele_Freqs_Initialize(Start_Freq[1:100], Mut_Drift_Dist, Start_Samp_Size[1:100])

# allele_freq_dist=sapply(indices, Bayes_Allele_Freq, Mut_Drift_Dist, .5, 26)

Allele_Freqs_Initialize = function(observed_freqs, prior_allele_freqs, allele_counts){
	nloci = length(observed_freqs)
	n_intervals = length(prior_allele_freqs)
	indices = c(1:n_intervals)
	start_frequencies = c()
	for(i in 1:nloci){
		obs_freq = observed_freqs[i]
		n_alleles = round(allele_counts[i])
		allele_freq_dist = sapply(indices, Bayes_Allele_Freq, prior_allele_freqs, obs_freq, n_alleles)
		allele_freq_dist = allele_freq_dist/sum(allele_freq_dist)
		c_allele_freq_dist = cumsum(allele_freq_dist)
		est_freq = (length(which(c_allele_freq_dist <= runif(1,0,1)))+1)/(n_intervals + 1)
		start_frequencies = c(start_frequencies, est_freq)
	}
return(start_frequencies)
}


	
#simulate Fisher-Wright
Fisher_Wright = function(ngen, init_freq, pop_size){
	n_loci = length(init_freq)
	start_sim = init_freq
	next_sim = start_sim
	for(iter in 1:ngen){
		for(loc in 1:n_loci){
			next_sim[loc]=rbinom(1,pop_size,start_sim[loc])/pop_size
		}
		start_sim = next_sim
	}
	return(next_sim)
}

#replicate nreps of Fisher-Wright model, return allele frequencies for each replicate
Reps_Fisher_Wright = function(nreps, ngen, init_freq, pop_size){
	freq_reps = c()
	for(i in 1:nreps){
		vout = Fisher_Wright(ngen, init_freq, pop_size)
		freq_reps = rbind(freq_reps, vout)
	}
	return(freq_reps)
}


Haploid_Sample = function(samp_size_dist, freq_dist){
	f_outputs = c()
	nsites = length(freq_dist)
	for(i in 1:nsites){
		f_outputs = c(f_outputs, rbinom(1, floor(samp_size_dist[i]), freq_dist[i])/samp_size_dist[i])
	}
	return(f_outputs)
}


Binom_Exclude_Fix = function(nsample, freq){
	x = 0
	while(x < 1 | x > (nsample-1)){
		x = rbinom(1, nsample, freq)
		#print(x)
	}
	return(x)
}

Haploid_Sample_No_Fix = function(samp_size_dist, freq_dist){
	f_outputs = c()
	nsites = length(freq_dist)
	for(i in 1:nsites){
		f_outputs = c(f_outputs, Binom_Exclude_Fix(floor(samp_size_dist[i]), freq_dist[i])/samp_size_dist[i])
	}
	return(f_outputs)
}


#mean of absolute frequency difference between initial and final frequencies
Avg_Abs_Freq_Diff = function(v1,v2){
	return(mean(abs(v1-v2)), na.rm=TRUE)
}


compare_allele_freq_change = function(v1, v2){
	return(mean(abs(v1 - v2),na.rm=TRUE))
}

#Fisher_Wright_With_End_Sampling(ngen, Start_Freq[1:100], popsize, Mut_Drift_Dist, Start_Samp_Size[1:100], End_Samp_Size[1:100], End_Raw_Samp_Size[1:100])
#Allele_Freqs_Initialize(Start_Freq[1:10], Mut_Drift_Dist, Start_Samp_Size[1:10])

#Bayesian approach - estimate initial frequency using equilibrium prior (Mut_Drift_Dist as prior)
#complete script includes sampling endpoints
Fisher_Wright_With_End_Sampling = function(n_gen, init_freq, pop_size, freq_prior, samp_size_init, chrom_size_final, samp_size_final){
	Start_Freq_New = Allele_Freqs_Initialize(init_freq, freq_prior, samp_size_init)
	Temp_Freq = Fisher_Wright(n_gen, Start_Freq_New, pop_size)
	Neff_Freq = Haploid_Sample(chrom_size_final, Temp_Freq)
	End_Freq = Haploid_Sample(samp_size_final, Neff_Freq)
	#End_Freq[is.nan(End_Freq)]=NA
	return(End_Freq)
}

#setwd("/raid10/mshpak/MuseumFlies/Counts/")
#Twentieth_2L = read.csv("AllCounts_Twentieth_Chr2L.txt",header=FALSE,sep='\t') 
#Twentieth_2R = read.csv("AllCounts_Twentieth_Chr2R.txt",header=FALSE,sep='\t')
Twentieth_3L = read.csv("AllCounts_Twentieth_Chr3L.txt",header=FALSE,sep='\t') 
#Twentieth_3R = read.csv("AllCounts_Twentieth_Chr3R.txt",header=FALSE,sep='\t')

#Twentieth_Autosome = rbind(Twentieth_2L, rbind(Twentieth_2R, rbind(Twentieth_3L, Twentieth_3R)))
Twentieth_Autosome = Twentieth_3L
n_alleles_1933 = rowSums(Twentieth_Autosome)
# at least 10 alleles - modify when used for 1800s
useful_sites_1933 = which(n_alleles_1933 > 9)

#for 3L - exclude telomere/centromere regions
useful_sites_1933 = useful_sites_1933[which(useful_sites_1933 > 600000 & useful_sites_1933 < 17700000)]
# others 2L: 0.5-17.5 Mb, 2R: 5.2-20.8 Mb, 3R: 6.9-26.6 Mb, 

find_max = function(v){
	tmax = max(v)
	return(which(v==tmax)[1])
}

#positions of most frequent allele in some population
max_pos = apply(Twentieth_Autosome, 1, find_max)
# this takes the max_pos position in each row, corresponding to the major allele count
max_allele_counts = Twentieth_Autosome[cbind(seq_len(nrow(Twentieth_Autosome)), max_pos)]

#identify polymorphic sites.
#take_sites = intersect(useful_sites_1933, which(max_allele_counts < n_alleles_1933 & (n_alleles_1933 - max_allele_counts > 1)))
# using < (n_alleles_1933 - 1) removes singletons from consideration
take_sites = intersect(useful_sites_1933, which(max_allele_counts <= (n_alleles_1933 - 1) & max_allele_counts >= ceiling(n_alleles_1933/2)))
#

allele_freqs_1933 = (max_allele_counts/n_alleles_1933)[take_sites]
samp_size_1933 = n_alleles_1933[take_sites]

# MODERN LUND SAMPLES

#Lund_2L = read.csv("AllCounts_ModernLund_Chr2L.txt",header=FALSE,sep='\t') 
#Lund_2R = read.csv("AllCounts_ModernLund_Chr2R.txt",header=FALSE,sep='\t')
Lund_3L = read.csv("AllCounts_ModernLund_Chr3L.txt",header=FALSE,sep='\t') 
#Lund_3R = read.csv("AllCounts_ModernLund_Chr3R.txt",header=FALSE,sep='\t')
#Lund_X = read.csv("AllCounts_ModernLund_ChrX.txt",header=FALSE,sep='\t')


# Lund raw allele counts - BEFORE poolseq effective sample size correction
#setwd("/home/mshpak/Lundflies/VCF/Pooled_Round1")
setwd("/raid10/mshpak/LundFlies/Count_Size/")
#Lund_RawCounts_2L = read.csv("Total_Counts_2L", header=FALSE, sep='\t')
#Lund_RawCounts_2R = read.csv("Total_Counts_2R", header=FALSE, sep='\t')
Lund_RawCounts_3L = read.csv("Total_Counts_3L", header=TRUE, sep=',')
#Lund_RawCounts_3L = read.csv("Total_Counts_3L", header=TRUE, sep=',')
#Lund_RawCounts_3R = read.csv("Total_Counts_3R", header=FALSE, sep='\t')
#Lund_RawCounts_X = read.csv("Total_Counts_X", header=FALSE, sep='\t')

Lund_RawCounts_3L = Lund_RawCounts_3L[,1:4]


# raw allele counts: actual size of sample before correction for effective sample size
#allele_counts_2L = rowSums(Lund_RawCounts_2L)
#allele_counts_2R = rowSums(Lund_RawCounts_2R)
allele_counts_3L = rowSums(Lund_RawCounts_3L)
#allele_counts_3R = rowSums(Lund_RawCounts_3R)
#allele_counts = c(allele_counts_2L, c(allele_counts_2R, c(allele_counts_3L, allele_counts_3R)))
allele_counts = allele_counts_3L

#Lund_Autosome = rbind(Lund_2L, rbind(Lund_2R, rbind(Lund_3L, Lund_3R)))
Lund_Autosome = Lund_3L
n_alleles_lund = rowSums(Lund_Autosome)
useful_sites_lund = which(n_alleles_lund > 19)

#lund allele counts at target loci based on polymorphic sites in 1933 sample
target_allele_count_lund = Lund_Autosome[cbind(seq_len(nrow(Lund_Autosome)),max_pos)]

allele_freq_lund = (target_allele_count_lund/n_alleles_lund)
allele_freq_lund[is.nan(allele_freq_lund)]=NA

#take_sites_new = intersect(take_sites, useful_sites_lund)

#allele frequencies at sites polymorphic in 1933
allele_freq_match = allele_freq_lund[take_sites]
samp_size_lund_match = n_alleles_lund[take_sites]
samp_size_lund_raw = allele_counts[take_sites]

#change_allele_freq = abs(allele_freq_match - allele_freqs_1933)

nreps = 100
ngens = 1230
#Start_Freq = allele_freqs_1933
#Final_Freq = allele_freq_match
#sample sizes per locus (read depth estimates)
#Start_Samp_Size = samp_size_1933
End_Samp_Size = samp_size_lund_match


#make sure read depth is nonzero at endpoint
# End_Raw_Samp_Size: read count, End_Samp_Size = effective sample size, i.e. effective read count
endpoint_sites = which(End_Samp_Size > 19)
Start_Freq = allele_freqs_1933[endpoint_sites]
Final_Freq = allele_freq_match[endpoint_sites]
Start_Samp_Size = samp_size_1933[endpoint_sites]
End_Samp_Size = samp_size_lund_match[endpoint_sites]
End_Raw_Samp_Size = samp_size_lund_raw[endpoint_sites]
nsites = length(End_Samp_Size)
#120 files, 240 autosomal chromosomes
End_Chromosomes = rep(240, nsites)

#test - exclude fixation 1/0
revised_endpoint = which(Final_Freq > 0 & Final_Freq < 1)
Start_Freq_Revised = Start_Freq[revised_endpoint]
End_Freq_Revised = Final_Freq[revised_endpoint]


#parallelize
library(foreach)
library(doParallel)
#numCores=64
numCores=20
registerDoParallel(numCores)

compare_allele_freq_change = function(v1, v2){
	return(mean(abs(v1 - v2),na.rm=TRUE))
}


compare_allele_freq_change(Start_Freq, Final_Freq)

#compare_allele_freq_change(Start_Freq, Final_Freq)
# for 3L 0.1652265

# Init_Freq = Allele_Freqs_Initialize(Start_Freq, Mut_Drift_Dist, Start_Samp_Size)
#Neff_Freq = Haploid_Sample(End_Samp_Size, Init_Freq)
#End_Freq = Haploid_Sample(End_Raw_Samp_Size, Neff_Freq)

popsize = 2000
nreps = 10
ngen = 1230

Mut_Drift_Dist = sapply(Freq_Interval, Mutation_Drift, dm_mut, popsize)
Mut_Drift_Dist = Mut_Drift_Dist/sum(Mut_Drift_Dist)

rubbish_2000 = foreach(i = 1:nreps, .combine=rbind) %dopar% {
	Fisher_Wright_With_End_Sampling(ngen, Start_Freq, popsize, Mut_Drift_Dist, Start_Samp_Size, End_Chromosomes, End_Samp_Size)
}

diffs_2000 = c()
for(i in 1:nreps){
	diffs_2000 = c(diffs_2000, compare_allele_freq_change(rubbish_2000[i,],Start_Freq))
}


write.csv(diffs_2000, file="chr3L_bayes_freq_differences_2K.csv", row.names=FALSE)
write.csv(rubbish_2000, file="chr3L_bayes_terminal_allele_freqs_2K.csv", row.names=FALSE)




popsize = 4000
nreps = 10
ngen = 1230
# try 820 for lower bound, i.e. 10 rather than 15 generations per year in Lund?

Mut_Drift_Dist = sapply(Freq_Interval, Mutation_Drift, dm_mut, popsize)
Mut_Drift_Dist = Mut_Drift_Dist/sum(Mut_Drift_Dist)

rubbish_4000 = foreach(i = 1:nreps, .combine=rbind) %dopar% {
	Fisher_Wright_With_End_Sampling(ngen, Start_Freq, popsize, Mut_Drift_Dist, Start_Samp_Size, End_Chromosomes, End_Samp_Size)
}

diffs_4000 = c()
for(i in 1:nreps){
	diffs_4000 = c(diffs_4000, compare_allele_freq_change(rubbish_4000[i,],Start_Freq))
}

write.csv(diffs_4000, file="chr3L_bayes_freq_differences_4K.csv", row.names=FALSE)
write.csv(rubbish_4000, file="chr3L_bayes_terminal_allele_freqs_4K.csv", row.names=FALSE)





popsize = 5000
nreps = 10
ngen = 1230
# try 820 for lower bound, i.e. 10 rather than 15 generations per year in Lund?

Mut_Drift_Dist = sapply(Freq_Interval, Mutation_Drift, dm_mut, popsize)
Mut_Drift_Dist = Mut_Drift_Dist/sum(Mut_Drift_Dist)


rubbish_5000 = foreach(i = 1:nreps, .combine=rbind) %dopar% {
	Fisher_Wright_With_End_Sampling(ngen, Start_Freq, popsize, Mut_Drift_Dist, Start_Samp_Size, End_Chromosomes, End_Samp_Size)
}

diffs_5000 = c()
for(i in 1:nreps){
	diffs_5000 = c(diffs_5000, compare_allele_freq_change(rubbish_5000[i,],Start_Freq))
}

write.csv(diffs_5000, file="chr3L_bayes_freq_differences_5K.csv", row.names=FALSE)
write.csv(rubbish_5000, file="chr3L_bayes_terminal_allele_freqs_5K.csv", row.names=FALSE)



popsize = 5400
nreps = 10
ngen = 1230

Mut_Drift_Dist = sapply(Freq_Interval, Mutation_Drift, dm_mut, popsize)
Mut_Drift_Dist = Mut_Drift_Dist/sum(Mut_Drift_Dist)

rubbish_5400 = foreach(i = 1:nreps, .combine=rbind) %dopar% {
	Fisher_Wright_With_End_Sampling(ngen, Start_Freq, popsize, Mut_Drift_Dist, Start_Samp_Size, End_Chromosomes, End_Samp_Size)
}

diffs_5400 = c()
for(i in 1:nreps){
	diffs_5400 = c(diffs_5400, compare_allele_freq_change(rubbish_5400[i,],Start_Freq))
}

write.csv(diffs_5400, file="chr3L_bayes_freq_differences_5400.csv", row.names=FALSE)
write.csv(rubbish_5400, file="chr3L_bayes_terminal_allele_freqs_5400.csv", row.names=FALSE)


popsize = 5600
nreps = 10
ngen = 1230

Mut_Drift_Dist = sapply(Freq_Interval, Mutation_Drift, dm_mut, popsize)
Mut_Drift_Dist = Mut_Drift_Dist/sum(Mut_Drift_Dist)

rubbish_5600 = foreach(i = 1:nreps, .combine=rbind) %dopar% {
	Fisher_Wright_With_End_Sampling(ngen, Start_Freq, popsize, Mut_Drift_Dist, Start_Samp_Size, End_Chromosomes, End_Samp_Size)
}

diffs_5600 = c()
for(i in 1:nreps){
	diffs_5600 = c(diffs_5600, compare_allele_freq_change(rubbish_5600[i,],Start_Freq))
}

write.csv(diffs_5600, file="chr3L_bayes_freq_differences_5600.csv", row.names=FALSE)
write.csv(rubbish_5600, file="chr3L_bayes_terminal_allele_freqs_5600.csv", row.names=FALSE)



popsize = 5800
nreps = 10
ngen = 1230

Mut_Drift_Dist = sapply(Freq_Interval, Mutation_Drift, dm_mut, popsize)
Mut_Drift_Dist = Mut_Drift_Dist/sum(Mut_Drift_Dist)

rubbish_5800 = foreach(i = 1:nreps, .combine=rbind) %dopar% {
	Fisher_Wright_With_End_Sampling(ngen, Start_Freq, popsize, Mut_Drift_Dist, Start_Samp_Size, End_Chromosomes, End_Samp_Size)
}

diffs_5800 = c()
for(i in 1:nreps){
	diffs_5800 = c(diffs_5800, compare_allele_freq_change(rubbish_5800[i,],Start_Freq))
}

write.csv(diffs_5800, file="chr3L_bayes_freq_differences_5800.csv", row.names=FALSE)
write.csv(rubbish_5800, file="chr3L_bayes_terminal_allele_freqs_5800.csv", row.names=FALSE)




popsize = 6200
nreps = 10
ngen = 1230

Mut_Drift_Dist = sapply(Freq_Interval, Mutation_Drift, dm_mut, popsize)
Mut_Drift_Dist = Mut_Drift_Dist/sum(Mut_Drift_Dist)

rubbish_6200 = foreach(i = 1:nreps, .combine=rbind) %dopar% {
	Fisher_Wright_With_End_Sampling(ngen, Start_Freq, popsize, Mut_Drift_Dist, Start_Samp_Size, End_Chromosomes, End_Samp_Size)
}

diffs_6200 = c()
for(i in 1:nreps){
	diffs_6200 = c(diffs_6200, compare_allele_freq_change(rubbish_6200[i,],Start_Freq))
}

write.csv(diffs_6200, file="chr3L_bayes_freq_differences_6200.csv", row.names=FALSE)
write.csv(rubbish_6200, file="chr3L_bayes_terminal_allele_freqs_6200.csv", row.names=FALSE)







