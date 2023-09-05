# simulate expected change in allele frequency under different effective population sizes
# 1809-1933 time interval

# redo for Chr 2R, 3L, exclude low recombination regions

#mutation rate per nucleotide (Keightley et al 2014)
dm_mut = 3*10**(-9)
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
#complete script includes sampling endpoints - note that for 1800s and 1933 museum flies, these are individual sequences, not poolseq
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
useful_sites_1933 = which(n_alleles_1933 > 9)

#for 3L
useful_sites_1933 = useful_sites_1933[which(useful_sites_1933 > 600000 & useful_sites_1933 < 17700000)]

#setwd("/raid10/mshpak/MuseumFlies/FAS1K/Fas_1800/Autosome/Lund_Only")
#Nineteenth_3L = read.csv("AllCounts_Nineteenth_LundOnly_Chr3L.txt",header=FALSE,sep='\t') 
Nineteenth_3L = read.csv("AllCounts_Nineteenth_LundOnly_Chr3L.txt",header=FALSE,sep='\t')
#Nineteenth_2R = read.csv("AllCounts_Nineteenth_Lund_Chr2R.txt",header=FALSE,sep='\t')

Nineteenth_Autosome = Nineteenth_3L
n_alleles_1800 = rowSums(Nineteenth_Autosome)
useful_sites_1800 = which(n_alleles_1800 > 3)

#for 3L
useful_sites_1800 = useful_sites_1800[which(useful_sites_1800 > 600000 & useful_sites_1800 < 17700000)]

find_max = function(v){
	tmax = max(v)
	return(which(v==tmax)[1])
}




#positions of most frequent allele in some population
max_pos = apply(Nineteenth_Autosome, 1, find_max)
# this takes the max_pos position in each row, corresponding to the major allele count
max_allele_counts = Nineteenth_Autosome[cbind(seq_len(nrow(Nineteenth_Autosome)), max_pos)]

#identify polymorphic sites.
#take_sites = intersect(useful_sites_1933, which(max_allele_counts < n_alleles_1933 & (n_alleles_1933 - max_allele_counts > 1)))
take_sites = intersect(useful_sites_1800, which(max_allele_counts <= (n_alleles_1800 - 1) & max_allele_counts >= ceiling(n_alleles_1800/2)))

allele_freqs_1800 = (max_allele_counts/n_alleles_1800)[take_sites]
samp_size_1800 = n_alleles_1800[take_sites]




#lund allele counts at target loci based on polymorphic sites in 1933 sample
target_allele_count_1933 = Twentieth_Autosome[cbind(seq_len(nrow(Twentieth_Autosome)),max_pos)]

allele_freq_1933 = (target_allele_count_1933/n_alleles_1933)
allele_freq_1933[is.nan(allele_freq_1933)]=NA

#take_sites_new = intersect(take_sites, useful_sites_lund)

#allele frequencies at sites polymorphic in 1933
allele_freq_match = allele_freq_1933[take_sites]
samp_size_match_1933 = n_alleles_1933[take_sites]

#change_allele_freq = abs(allele_freq_match - allele_freqs_1933)

nreps = 10
ngens = 1860
#Start_Freq = allele_freqs_1933
#Final_Freq = allele_freq_match
#sample sizes per locus (read depth estimates)
#Start_Samp_Size = samp_size_1933
End_Samp_Size = samp_size_match_1933


#make sure read depth is nonzero at endpoint
endpoint_sites = which(End_Samp_Size > 9)
Start_Freq = allele_freqs_1800[endpoint_sites]
Final_Freq = allele_freq_match[endpoint_sites]
Start_Samp_Size = samp_size_1800[endpoint_sites]
End_Samp_Size = samp_size_match_1933[endpoint_sites]
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
# for 3L 0.2431251

#Allele_Freqs_Initialize = function(observed_freqs, prior_allele_freqs, allele_counts){

popsize = 4000
nreps = 10
ngen = 1860

Mut_Drift_Dist = sapply(Freq_Interval, Mutation_Drift, dm_mut, popsize)
Mut_Drift_Dist = Mut_Drift_Dist/sum(Mut_Drift_Dist)


rubbish_4000 = foreach(i = 1:nreps, .combine=rbind) %dopar% {
	Fisher_Wright_With_End_Sampling(ngen, Start_Freq, popsize, Mut_Drift_Dist, Start_Samp_Size, End_Chromosomes, End_Samp_Size)
}

diffs_4000 = c()
for(i in 1:nreps){
	diffs_4000 = c(diffs_4000, compare_allele_freq_change(rubbish_4000[i,],Start_Freq))
}

write.csv(diffs_4000, file="chr3L_bayes_freq_differences_1800_4K.csv", row.names=FALSE)
write.csv(rubbish_4000, file="chr3L_bayes_terminal_allele_freqs_1800_4K.csv", row.names=FALSE)



popsize = 3000
nreps = 10
ngen = 1860

Mut_Drift_Dist = sapply(Freq_Interval, Mutation_Drift, dm_mut, popsize)
Mut_Drift_Dist = Mut_Drift_Dist/sum(Mut_Drift_Dist)


rubbish_3000 = foreach(i = 1:nreps, .combine=rbind) %dopar% {
	Fisher_Wright_With_End_Sampling(ngen, Start_Freq, popsize, Mut_Drift_Dist, Start_Samp_Size, End_Chromosomes, End_Samp_Size)
}

diffs_3000 = c()
for(i in 1:nreps){
	diffs_3000 = c(diffs_3000, compare_allele_freq_change(rubbish_3000[i,],Start_Freq))
}

write.csv(diffs_3000, file="chr3L_bayes_freq_differences_1800_3K.csv", row.names=FALSE)
write.csv(rubbish_3000, file="chr3L_bayes_terminal_allele_freqs_1800_3K.csv", row.names=FALSE)



popsize = 3400
nreps = 10
ngen = 1860
# try 820 for lower bound, i.e. 10 rather than 15 generations per year in Lund?

Mut_Drift_Dist = sapply(Freq_Interval, Mutation_Drift, dm_mut, popsize)
Mut_Drift_Dist = Mut_Drift_Dist/sum(Mut_Drift_Dist)


rubbish_3400 = foreach(i = 1:nreps, .combine=rbind) %dopar% {
	Fisher_Wright_With_End_Sampling(ngen, Start_Freq, popsize, Mut_Drift_Dist, Start_Samp_Size, End_Chromosomes, End_Samp_Size)
}

diffs_3400 = c()
for(i in 1:nreps){
	diffs_3400 = c(diffs_3400, compare_allele_freq_change(rubbish_3400[i,],Start_Freq))
}

write.csv(diffs_3400, file="chr3L_bayes_freq_differences_1800_3400.csv", row.names=FALSE)
write.csv(rubbish_3400, file="chr3L_bayes_terminal_allele_freqs_1800_3400.csv", row.names=FALSE)


popsize = 3600
nreps = 10
ngen = 1860
# try 820 for lower bound, i.e. 10 rather than 15 generations per year in Lund?

Mut_Drift_Dist = sapply(Freq_Interval, Mutation_Drift, dm_mut, popsize)
Mut_Drift_Dist = Mut_Drift_Dist/sum(Mut_Drift_Dist)


rubbish_3600 = foreach(i = 1:nreps, .combine=rbind) %dopar% {
	Fisher_Wright_With_End_Sampling(ngen, Start_Freq, popsize, Mut_Drift_Dist, Start_Samp_Size, End_Chromosomes, End_Samp_Size)
}

diffs_3600 = c()
for(i in 1:nreps){
	diffs_3600 = c(diffs_3600, compare_allele_freq_change(rubbish_3600[i,],Start_Freq))
}

write.csv(diffs_3600, file="chr3L_bayes_freq_differences_1800_3600.csv", row.names=FALSE)
write.csv(rubbish_3600, file="chr3L_bayes_terminal_allele_freqs_1800_3600.csv", row.names=FALSE)


popsize = 3800
nreps = 10
ngen = 1860
# try 820 for lower bound, i.e. 10 rather than 15 generations per year in Lund?

Mut_Drift_Dist = sapply(Freq_Interval, Mutation_Drift, dm_mut, popsize)
Mut_Drift_Dist = Mut_Drift_Dist/sum(Mut_Drift_Dist)


rubbish_3800 = foreach(i = 1:nreps, .combine=rbind) %dopar% {
	Fisher_Wright_With_End_Sampling(ngen, Start_Freq, popsize, Mut_Drift_Dist, Start_Samp_Size, End_Chromosomes, End_Samp_Size)
}

diffs_3800 = c()
for(i in 1:nreps){
	diffs_3800 = c(diffs_3800, compare_allele_freq_change(rubbish_3800[i,],Start_Freq))
}

write.csv(diffs_3800, file="chr3L_bayes_freq_differences_1800_3800.csv", row.names=FALSE)
write.csv(rubbish_3800, file="chr3L_bayes_terminal_allele_freqs_1800_3800.csv", row.names=FALSE)



