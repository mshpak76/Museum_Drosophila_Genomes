# redo IBD calculations for all pairs
# consider pairs H4, H6, H10, H11, H12, H13 for 1800s, all for 1933
# i.e. H1, H3, H7  and 10-24

setwd("/raid10/mshpak/MuseumFlies/IBD_Calculations")

chrom = "Chr3R" #change as needed
outfile_name = paste("Full_Mat_", chrom, sep = '')
outfile_name = paste(outfile_name, ".csv", sep = '')

#write.csv(Full_Mat, file = outfile_name, row.names=FALSE)
#
Full_Mat = read.csv(outfile_name, header=TRUE)
Full_Mat = Full_Mat[,1:24]

#chr_window_filename = paste("windows_ZI10000_", chrom, sep = '')
chr_window_filename = paste("windows_ZI25000_", chrom, sep = '')
chr_window_filename = paste(chr_window_filename, ".txt", sep = '')



#with IBD:
# H4, H6, H12, H13
# H1, H3, H7, H8, H14, H15, H16, H18, H23
# small overall H5
#background_exclude = c(1,3,4,7,10,11,12,13,14,15,16,18,23)


#3R
elements_1800 = c(3,1,4,7,6,5)
pairs_1800 = c()
for(i in 1:(length(elements_1800)-1)){
	for(j in (i+1):length(elements_1800)){
		pairs_1800 = rbind(pairs_1800, c(elements_1800[i],elements_1800[j]))
	}
}


# 1933 priority masking order H3, 7, 1, 16, 15, 21, 8, 22, 18, 14, 17, 19, 20, 24, 23
elements_1933 = c(10, 11, 12, 16, 15, 21, 13, 22, 18, 14, 17, 14, 20, 24, 23) 
pairs_1933 = c()
for(i in 1:(length(elements_1933)-1)){
	for(j in (i+1):length(elements_1933)){
		pairs_1933 = rbind(pairs_1933, c(elements_1933[i], elements_1933[j]))
	}
}

# when calculating average distance between an individual chromosome and others in the population
# for comparison to outgroup, we exclude those individuals with obvious pairwise IBD based on
# distance metrics so as not to deflate genetic distance

assoc_exclude_1800 = list()
assoc_exclude_1800[[1]] = c(1,3,4)
assoc_exclude_1800[[3]] = c(3,1,4,7)
assoc_exclude_1800[[4]] = c(4,1,3,5,7)
assoc_exclude_1800[[5]] = c(5,4,7)
assoc_exclude_1800[[6]] = c(6)
assoc_exclude_1800[[7]] = c(7,3,4,5)

assoc_exclude_1933 = list()
assoc_exclude_1933[[10]] = c(10,11,12,19,24)
assoc_exclude_1933[[11]] = c(11,10,19,24)
assoc_exclude_1933[[12]] = c(12,10)
assoc_exclude_1933[[13]] = c(8)
assoc_exclude_1933[[14]] = c(14,17,18,19)
assoc_exclude_1933[[15]] = c(15,20)
assoc_exclude_1933[[16]] = c(16,23)
assoc_exclude_1933[[17]] = c(17,14)
assoc_exclude_1933[[18]] = c(18,14,19)
assoc_exclude_1933[[19]] = c(19,10,11,14,18,22)
assoc_exclude_1933[[20]] = c(20,15,23)
assoc_exclude_1933[[21]] = c(21,24)
assoc_exclude_1933[[22]] = c(22,19)
assoc_exclude_1933[[23]] = c(23,16,20)
assoc_exclude_1933[[24]] = c(24,10,11,21)


#try alternative seq prioritizing - consider order of low quality
seq_to_mask_prioritize = function(mat, low_qual_samples){
	temp_mat = mat
	to_mask = c()
	all_instances = c(temp_mat[,1],temp_mat[,2])
	if(any(low_qual_samples %in% all_instances)){
		iter = 0
		lq = low_qual_samples[which(low_qual_samples %in% all_instances)]
		while(length(temp_mat[,1] > 0) & any(lq %in% all_instances)){
			lq_samp = lq[1]
			low_qual_pos = which(lq_samp %in% temp_mat[,1] | lq_samp %in% temp_mat[,2])
			low_qual_add = lq[low_qual_pos]
			#to_mask = c(to_mask,low_qual_add)
			posits = which(temp_mat[,1] %in% low_qual_add | temp_mat[,2] %in% low_qual_add)
			if(length(temp_mat[,1]) > 1){
				to_mask = c(to_mask, lq_samp)
			}else{
				to_mask = c(to_mask, temp_mat[1,1])
			}
			lq = lq[2:length(lq)]
			temp_mat = temp_mat[-posits,,drop=FALSE]
			all_instances = c(temp_mat[,1],temp_mat[,2])
			iter = iter + 1
		}		
	}
	iter = 0
	while(length(temp_mat[,1]>0)){
		all_instances = c(temp_mat[,1],temp_mat[,2])
		#number of times a sample appears
		occurrences = sort(table(all_instances))
		max_occur = max(occurrences)
		which_max_occur = which(occurrences == max_occur)
		max_instances = as.numeric(names(which_max_occur)[1])
		posits = which(temp_mat[,1]==max_instances | temp_mat[,2]==max_instances)
		if(length(temp_mat[,1]) > 1){
			to_mask = c(to_mask, max_instances)
		}else{
			to_mask = c(to_mask, temp_mat[1,1])
		}
		temp_mat = temp_mat[-posits,,drop=FALSE]
		iter = iter + 1
	}
	return(to_mask)
}


ibd_prioritize_1800 = seq_to_mask_prioritize(ibd_pairs_1800,c())
ibd_prioritize_1933 = seq_to_mask_prioritize(ibd_pairs_1933,c())

#keep
for(i in assoc_exclude_1800){
	print(museum_samples[i])
}

for(i in assoc_exclude_1933){
	print(museum_samples[i])
}

seq_ham_dist = function(v1,v2){
	return(sum(v1 != v2, na.rm = TRUE))
}

problem_pairs_1800 = pairs_1800
problem_pairs_1933 = pairs_1933

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

for(i in 1:length(problem_pairs_1800[,1])){
	print(c(museum_samples[problem_pairs_1800[i,1]],museum_samples[problem_pairs_1800[i,2]]))
}
 
window_intervals = read.csv(chr_window_filename, sep = '\t', header=FALSE)

dict_1800 = list()
for(i in 1:length(problem_pairs_1800[,1])){
	dict_1800[[museum_samples[problem_pairs_1800[i,1]]]]=c()
}

dict_1933 = list()
for(i in 1:length(problem_pairs_1933[,1])){
	dict_1933[[museum_samples[problem_pairs_1933[i,1]]]]=c()
}


#low_qual_1800 = c(3,1,4)
#low_qual_1933 = c(10,11)

#low_qual_1800 = c(low_qual_1800, setdiff(ibd_prioritize_1800, low_qual_1800))
#low_qual_1933 = c(low_qual_1933, setdiff(ibd_prioritize_1933, low_qual_1933))

#low_qual_1800 = c(3,1,4,5,6,7)
# use this solely for case where there's strict priority based on coverage, including high-quality sequences
#low_qual_1933 = c(10,11,12,15,16,14,13,21,22,18,19,20,17,24,23)


low_qual_1800 = elements_1800
low_qual_1933 = elements_1933

all_1800 = unique(c(problem_pairs_1800[,1],problem_pairs_1800[,2]))
all_1933 = unique(c(problem_pairs_1933[,1],problem_pairs_1933[,2]))


# this one prioritizes masking of low quality sequences, those with the greatest amount of masking, and highest occurrence, in that order

#create correction factors based on average distances to outside populations
window_correction_factors = c() #per window per sample
window_pi_1800 = c() # per window across all samples
window_pi_1933 = c()
for(i in 1:length(window_intervals[,1])){
	start = window_intervals[i,1]+1
	end = window_intervals[i,2]+1
	part_mat = Full_Mat[start:end,]
	Dist_Part_Mat = comp_mat_distances(part_mat,museum_samples)
	correct_factors = rep(1,length(museum_samples))
	avg_win_dist_1800 = c()
	for(j in all_1800){
		vtemp = Dist_Part_Mat[j, all_1800]
		vtemp = vtemp[-assoc_exclude_1800[[j]]]
		med_dist = median(vtemp, na.rm = TRUE)
		correct_factors[j] = med_dist
		avg_win_dist_1800 = c(avg_win_dist_1800, med_dist)
	}
	pop_avg_1800 = median(avg_win_dist_1800, na.rm = TRUE)
	window_pi_1800 = c(window_pi_1800, pop_avg_1800)
	correct_factors[all_1800] = correct_factors[all_1800]/pop_avg_1800
	avg_win_dist_1933 = c()
	for(j in all_1933){
		vtemp = Dist_Part_Mat[j, all_1933]
		vtemp = vtemp[-assoc_exclude_1933[[j]]]
		med_dist = median(vtemp, na.rm = TRUE)
		correct_factors[j] = med_dist
		avg_win_dist_1933 = c(avg_win_dist_1933, med_dist)
	}
	pop_avg_1933 = median(avg_win_dist_1933, na.rm = TRUE)
	window_pi_1933 = c(window_pi_1933, pop_avg_1933)
	correct_factors[all_1933] = correct_factors[all_1933]/pop_avg_1933
	window_correction_factors = rbind(window_correction_factors, correct_factors)
}


for(i in 1:length(window_intervals[,1])){
#for(i in 10:11){
	start = window_intervals[i,1]+1
	end = window_intervals[i,2]+1
	part_mat = Full_Mat[start:end,]
	Dist_Part_Mat = comp_mat_distances(part_mat,museum_samples) 
	temp_pairs_1800 = c()
	for(j in 1:length(problem_pairs_1800[,1])){
		first_samp = problem_pairs_1800[j,1]
		second_samp = problem_pairs_1800[j,2]
		dist_pair = Dist_Part_Mat[first_samp, second_samp]
		if(dist_pair < 0.875*window_pi_1800[i]*window_correction_factors[i,first_samp]*window_correction_factors[i,second_samp]){
			temp_pairs_1800 = rbind(temp_pairs_1800, problem_pairs_1800[j,])
		}
	}
	seq_mask_1800 = seq_to_mask_prioritize(temp_pairs_1800, low_qual_1800)
	for(val in seq_mask_1800){
		dict_1800[[museum_samples[val]]] = unique(c(dict_1800[[museum_samples[val]]], i))
	}
	## 1933 samples
	temp_pairs_1933 = c()
	for(j in 1:length(problem_pairs_1933[,1])){
		first_samp = problem_pairs_1933[j,1]
		second_samp = problem_pairs_1933[j,2]
		dist_pair = Dist_Part_Mat[first_samp, second_samp]
		if(dist_pair < 0.875*window_pi_1933[i]*window_correction_factors[i,first_samp]*window_correction_factors[i,second_samp]){
			temp_pairs_1933 = rbind(temp_pairs_1933, problem_pairs_1933[j,])
		}
	}
	seq_mask_1933 = seq_to_mask_prioritize(temp_pairs_1933, low_qual_1933)
	for(val in seq_mask_1933){
		dict_1933[[museum_samples[val]]] = unique(c(dict_1933[[museum_samples[val]]], i))
		
	}
}

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
write.csv(To_Remove, file="Remove_Windows_3R.csv", row.names=FALSE)


to_remove_s = read.csv("Remove_Windows_3R.csv")
for(i in 1:length(to_remove_s[,1])){
	samp_s = to_remove_s[i,1]
	windows_s = to_remove_s[i,]
	windows_s = windows_s[-1]
	windows_s = windows_s[which(windows_s != "NA")]
	print(samp_s)
	print(length(windows_s))
}
