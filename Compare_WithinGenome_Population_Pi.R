# identify possible IBD regions within genomes
# compare heterozygosity in a region to median pi in background

#for autosomes
#1800s: H10, H11, H13, H9, H25
#1933: H1, H8, H14-24

setwd("/raid10/mshpak/MuseumFlies/IBD_Calculations")

chrom = "Chr3L" #change as needed
outfile_name = paste("Full_Mat_", chrom, sep = '')
outfile_name = paste(outfile_name, ".csv", sep = '')


Full_Mat = read.csv(outfile_name, header=TRUE)
Full_Mat = Full_Mat[,1:24]

#chr_window_filename = paste("windows_ZI10000_", chrom, sep = '')
chr_window_filename = paste("windows_ZI25000_", chrom, sep = '')
chr_window_filename = paste(chr_window_filename, ".txt", sep = '')

window_intervals = read.csv(chr_window_filename, sep = '\t', header=FALSE)

#3L
elements_1800 = c(3,1,4,6,7,5)
pairs_1800 = c()
for(i in 1:(length(elements_1800)-1)){
	for(j in (i+1):length(elements_1800)){
		pairs_1800 = rbind(pairs_1800, c(elements_1800[i],elements_1800[j]))
	}
}


# 1933 priority masking order H3, 7, 1, 16, 15, 21, 8, 22, 14, 17, 18, 19, 20, 24, 23
elements_1933 = c(10, 11, 12, 16, 15, 21, 13, 22, 14, 17, 18, 19, 20, 24, 23)
pairs_1933 = c()
for(i in 1:(length(elements_1933)-1)){
	for(j in (i+1):length(elements_1933)){
		pairs_1933 = rbind(pairs_1933, c(elements_1933[i], elements_1933[j]))
	}
}

assoc_exclude_1800 = list()
assoc_exclude_1800[[1]] = c(1,3,4)
assoc_exclude_1800[[3]] = c(3,1,4,5,7)
assoc_exclude_1800[[4]] = c(4,1,3,6)
assoc_exclude_1800[[5]] = c(5,3,4)
assoc_exclude_1800[[6]] = c(6,4)
assoc_exclude_1800[[7]] = c(7,3)
assoc_exclude_1800[[8]] = c(8,9)
assoc_exclude_1800[[9]] = c(9,8)


assoc_exclude_1933 = list()
assoc_exclude_1933[[10]] = c(10,12,19)
assoc_exclude_1933[[11]] = c(11,24)
assoc_exclude_1933[[12]] = c(12,10,13,14)
assoc_exclude_1933[[13]] = c(13,14)
assoc_exclude_1933[[14]] = c(14,13,22)
assoc_exclude_1933[[15]] = c(15,23)
assoc_exclude_1933[[16]] = c(16,23)
assoc_exclude_1933[[17]] = c(17)
assoc_exclude_1933[[18]] = c(18,19)
assoc_exclude_1933[[19]] = c(19,10,18)
assoc_exclude_1933[[20]] = c(20)
assoc_exclude_1933[[21]] = c(21)
assoc_exclude_1933[[22]] = c(22,13)
assoc_exclude_1933[[23]] = c(23,15,16)
assoc_exclude_1933[[24]] = c(24,11)



seq_ham_dist = function(v1,v2){
	return(sum(v1 != v2, na.rm = TRUE))
}


# average distance in relation to non-na sites among pairs
pairwise_pi = function(v1,v2){
	nsites = which(is.na(v1)==FALSE & is.na(v2)==FALSE)
	sdist = seq_ham_dist(v1,v2)
	return(sdist/length(nsites))
}


comp_mat_distances = function(M1){
	Dmat = c()
	for(i in 1:length(M1[1,])){
		vtemp = c()
		for(j in 1:length(M1[1,])){
			dist = pairwise_pi(M1[,i],M1[,j])
			vtemp = c(vtemp, dist)
		}
		Dmat = rbind(Dmat, vtemp)
	}
	return(Dmat)
}

comp_mat_distances = function(M1,samples){
	Dmat = c()
	for(i in 1:length(samples)){
		vtemp = c()
		for(j in 1:length(M1[1,])){
			dist = pairwise_pi(M1[,i],M1[,j])
			vtemp = c(vtemp, dist)
		}
		Dmat = rbind(Dmat, vtemp)
	}
	return(Dmat)
}


Hap1800 = c('H4','H5','H6','H12')
Hap1933 = c('H3','H7')
Dip1800 = c('H10','H11','H13','H9','H25')
Dip1933 = c('H1','H8','H14','H15','H16','H17','H18','H19','H20','H21','H22','H23','H24')
museum_samples = c(Hap1800,Dip1800,Hap1933,Dip1933)

all_1800 = unique(elements_1800)
all_1933 = unique(elements_1933)


#window_correction_factors = c() #per window per sample
window_pi_1800 = c() # per window across all samples
window_pi_1933 = c()
for(i in 1:length(window_intervals[,1])){
	start = window_intervals[i,1]+1
	end = window_intervals[i,2]+1
	part_mat = Full_Mat[start:end,]
	Dist_Part_Mat = comp_mat_distances(part_mat,museum_samples)
	#correct_factors = rep(1,length(museum_samples))
	avg_win_dist_1800 = c()
	for(j in all_1800){
		vtemp = Dist_Part_Mat[j, all_1800]
		vtemp = vtemp[-assoc_exclude_1800[[j]]]
		med_dist = median(vtemp, na.rm = TRUE)
		#correct_factors[j] = med_dist
		avg_win_dist_1800 = c(avg_win_dist_1800, med_dist)
	}
	pop_avg_1800 = median(avg_win_dist_1800, na.rm = TRUE)
	window_pi_1800 = c(window_pi_1800, pop_avg_1800)
	#correct_factors[all_1800] = correct_factors[all_1800]/pop_avg_1800
	avg_win_dist_1933 = c()
	for(j in all_1933){
		vtemp = Dist_Part_Mat[j, all_1933]
		vtemp = vtemp[-assoc_exclude_1933[[j]]]
		med_dist = median(vtemp, na.rm = TRUE)
		#correct_factors[j] = med_dist
		avg_win_dist_1933 = c(avg_win_dist_1933, med_dist)
	}
	pop_avg_1933 = median(avg_win_dist_1933, na.rm = TRUE)
	window_pi_1933 = c(window_pi_1933, pop_avg_1933)
	#correct_factors[all_1933] = correct_factors[all_1933]/pop_avg_1933
	#window_correction_factors = rbind(window_correction_factors, correct_factors)
}


#number of distinct nucleotides (0,1 or 2)
nuc_count = function(v){
	return(length(which(v>0)))
}

heterozygosity = function(M1){
	vsum1 = rowSums(M1)
	to_keep = which(vsum1 > 0)
	mat_temp = M1[to_keep,]
	allele_counts = apply(mat_temp, 1, nuc_count)
	het = length(which(allele_counts > 1))/length(to_keep)
	return(het)
}
	

setwd("/raid10/mshpak/MuseumFlies/IndivCounts/MuseumIndivs")

#Count_File = read.csv(indiv_count_name, sep='\t', header=FALSE)


dict_1800 = list()
dict_1933 = list()

samp_1933 = c("H1","H8","H14","H15","H16","H17","H18","H19","H20","H21","H22","H23","H24")
samp_1800 = c("H10","H11","H13","H9","H25")
#samp_1933 = c("H1","H8")
for(samp in samp_1933){
	indiv_count_name = paste("IndivCounts_", samp, sep='')
	indiv_count_name = paste(indiv_count_name, "_NoInv_", sep='')
	indiv_count_name = paste(indiv_count_name, chrom, sep='')
	indiv_count_name = paste(indiv_count_name, ".txt", sep='')
	Count_File = read.csv(indiv_count_name, sep='\t', header=FALSE)
	het_windows = c()
	to_mask = c()
	for(i in 1:length(window_intervals[,1])){
		start = window_intervals[i,1]+1
		end = window_intervals[i,2]+1
		Part_Mat = Count_File[start:end,]
		het_temp = heterozygosity(Part_Mat)
		het_windows = c(het_windows, het_temp)
		if(het_temp < window_pi_1933/2){
			to_mask = c(to_mask, i)
		}
	}
	dict_1933[[samp]] = to_mask
}

for(samp in samp_1800){
	indiv_count_name = paste("IndivCounts_", samp, sep='')
	indiv_count_name = paste(indiv_count_name, "_NoInv_", sep='')
	indiv_count_name = paste(indiv_count_name, chrom, sep='')
	indiv_count_name = paste(indiv_count_name, ".txt", sep='')
	Count_File = read.csv(indiv_count_name, sep='\t', header=FALSE)
	het_windows = c()
	to_mask = c()
	for(i in 1:length(window_intervals[,1])){
		start = window_intervals[i,1]+1
		end = window_intervals[i,2]+1
		Part_Mat = Count_File[start:end,]
		het_temp = heterozygosity(Part_Mat)
		het_windows = c(het_windows, het_temp)
		if(het_temp < window_pi_1800/2){
			to_mask = c(to_mask, i)
		}
	}
	dict_1800[[samp]] = to_mask
}

#


max_len = length(window_intervals[,1])

To_Remove_1800 = c()
for(samp in names(dict_1800)){
	tempv = c(samp, dict_1800[[samp]])
	length(tempv)=max_len
	To_Remove_1800 = rbind(To_Remove_1800, tempv)
}

To_Remove_1933 = c()
for(samp in names(dict_1933)){
	tempv = c(samp, dict_1933[[samp]])
	length(tempv)=max_len
	To_Remove_1933 = rbind(To_Remove_1933, tempv)
}

To_Remove = rbind(To_Remove_1800, To_Remove_1933)


# change this to appropriate chromosome
write.csv(To_Remove, file="WithinGenome_Remove_Windows_3L.csv", row.names=FALSE)
#
	
