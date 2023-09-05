setwd("/raid10/mshpak/MuseumFlies/IBD_Calculations")

#takes output tables from Identify_IBD regions scripts and masks the corresponding coordinates from
#individual sample/chromosome count files

chrom = "2L"

win_remove_filename = paste("Remove_Windows_", chrom, sep='')
#win_remove_filename = paste("Passau_Remove_Windows_", chrom, sep='')
win_remove_filename = paste(win_remove_filename, ".csv", sep='')

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

#count file
setwd("/raid10/mshpak/MuseumFlies/IndivCounts/MuseumIndivs")

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

samples = names(coord_list)
#samples = samples[1]
for(samp in samples){
	count_file_name = paste("IndivCounts_", samp, sep ='')
	count_file_name = paste(count_file_name, "_NoInv_Chr", sep = '')
	count_file_name = paste(count_file_name, chrom, sep = '')
	count_file_name = paste(count_file_name, ".txt", sep = '')
	Count_Mat = read.csv(count_file_name, header=FALSE, sep='\t')
	To_Remove = coord_list[[samp]]
	for(i in 1:length(To_Remove[,1])){
		start_pos = To_Remove[i,1]
		end_pos = To_Remove[i,2]
		nsize = end_pos - start_pos + 1
		temp_mat = matrix(0, nsize, 4)
		Count_Mat[start_pos:end_pos, ] = temp_mat
	}
	outname = paste("New_IndivCounts_", samp, sep = '')
	outname = paste(outname, "_Chr", sep = '')
	outname = paste(outname, chrom, sep = '')
	outname = paste(outname, ".txt", sep = '')
	write.table(Count_Mat, outname, row.names=FALSE, col.names=FALSE, sep = '\t')
}
		
#			
		

setwd("/raid10/mshpak/MuseumFlies/IBD_Calculations")
