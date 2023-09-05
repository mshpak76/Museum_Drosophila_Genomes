setwd("/raid10/mshpak/MuseumFlies/IBD_Calculations")

#mask homozygous regions in the count files

chrom = "X"

Windows_To_Remove = read.csv(win_remove_filename, header=TRUE)

#window_intervals_filename = paste("windows_ZI10000_Chr", chrom, sep='')
window_intervals_filename = paste("windows_ZI25000_Chr", chrom, sep='')
window_intervals_filename = paste(window_intervals_filename, ".txt", sep='')

Windows_Coords = read.csv(window_intervals_filename, sep= '\t', header=FALSE)
total_windows = length(Windows_Coords[,1])

#setwd("/raid10/mshpak/MuseumFlies/IBD_Calculations/Windows5000")
#half_window_intervals_filename = paste("windows_ZI5000_Chr", chrom, sep='')
setwd("/raid10/mshpak/MuseumFlies/IBD_Calculations/Windows12500")
half_window_intervals_filename = paste("windows_ZI12500_Chr", chrom, sep='')
half_window_intervals_filename = paste(half_window_intervals_filename, ".txt", sep='')

Half_Windows_Coords = read.csv(half_window_intervals_filename, sep= '\t', header=FALSE)

setwd("/raid10/mshpak/MuseumFlies/IndivCounts/MuseumIndivs")

win_remove_filename = paste("WithinGenome_Remove_Windows_", chrom, sep='')
win_remove_filename = paste(win_remove_filename, ".csv", sep='')
Windows_To_Remove = read.csv(win_remove_filename, header=TRUE)

coord_list = list()
for(i in 1:length(Windows_To_Remove[,1])){
	sample = toString(Windows_To_Remove[i,1])
	coord_list[[sample]] = c()
	no_na_elements = length(which(Windows_To_Remove[i,] != "NA"))
	for(j in 2:no_na_elements){
		get_range = as.numeric(Windows_To_Remove[[i,j]])
		start_end = as.numeric(Windows_Coords[get_range,]) + 1
		#sample = toString(Windows_To_Remove[i,1])
		coord_list[[sample]] = rbind(coord_list[[sample]],start_end)
		#print(sample)
		if(get_range > 1){
			# half of previous interval
			#get_range = as.numeric(Windows_To_Remove[[i,j]])
			get_range_half = 2*(get_range - 1)
			start_end = as.numeric(Half_Windows_Coords[get_range_half,]) + 1
			coord_list[[sample]] = rbind(coord_list[[sample]], start_end)
		}
		if(get_range < total_windows){
			# half of next interval
			#get_range = as.numeric(Windows_To_Remove[[i,j]])
			get_range_half = 2*get_range + 1
			start_end = as.numeric(Half_Windows_Coords[get_range_half,]) + 1
			coord_list[[sample]] = rbind(coord_list[[sample]], start_end)
		}
	}
}

make_homozygous = function(v){
	vout = v
	rsum = sum(v,na.rm=TRUE)
	if(rsum > 0){
		alleles = which(v > 0)
		if(length(alleles) > 1){
			vout[sample(alleles,1)] = 0
		}else{
			vout[alleles[1]] = 1
		}
	}
	return(vout)
}

mask_homozygous = function(M){
	NewMat = t(apply(M, 1, make_homozygous))
	return(NewMat)
}			
	

samples = names(coord_list)
#samples = samples[1]
for(samp in samples){
	count_file_name = paste("New_IndivCounts_", samp, sep ='')
	count_file_name = paste(count_file_name, "_Chr", sep = '')
	count_file_name = paste(count_file_name, chrom, sep = '')
	count_file_name = paste(count_file_name, ".txt", sep = '')
	Count_Mat = read.csv(count_file_name, header=FALSE, sep='\t')
	To_Remove = coord_list[[samp]]
	for(i in 1:length(To_Remove[,1])){
		start_pos = To_Remove[i,1]
		end_pos = To_Remove[i,2]
		Submat = Count_Mat[start_pos:end_pos,]
		To_Substitute = mask_homozygous(Submat)
		Count_Mat[start_pos:end_pos, ] = To_Substitute
	}
	outname = paste("IBD_New_IndivCounts_", samp, sep = '')
	outname = paste(outname, "_Chr", sep = '')
	outname = paste(outname, chrom, sep = '')
	outname = paste(outname, ".txt", sep = '')
	write.table(Count_Mat, outname, row.names=FALSE, col.names=FALSE, sep = '\t')
}
