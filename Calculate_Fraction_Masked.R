setwd("/raid10/mshpak/MuseumFlies/IndivCounts/MuseumIndivs")

samples = c("H1","H3","H4","H6","H7","H8","H9","H10","H11","H12","H13","H14","H15","H16","H17","H18","H19","H20","H21","H22","H23","H24")
chroms = c("Chr2L", "Chr2R","Chr3L", "Chr3R", "ChrX")
#chroms = c("ChrX","Chr3R", "Chr3L")
#chroms = c("ChrX")

summary_mat = c()
for(chrom in chroms){
	for(samp in samples){
		oldname = paste("IndivCounts_", samp, sep='')
		oldname = paste(oldname, "_NoInv_", sep='')
		oldname = paste(oldname, chrom, sep = '')
		oldname = paste(oldname, ".txt", sep = '')
		newname = paste("New_IndivCounts_", samp, sep = '')
		newname = paste(newname, "_", sep = '')
		newname = paste(newname, chrom, sep = '')
		newname = paste(newname, ".txt", sep = '')
		Old = read.csv(oldname, sep='\t', header=FALSE)
		New = read.csv(newname, sep='\t', header=FALSE)
		rsum_old = rowSums(Old)
		rsum_new = rowSums(New)
		masked = 1-length(which(rsum_new > 0))/length(which(rsum_old > 0))
		print(chrom)
		print(samp)
		print(masked)
		dat_temp = c(samp, chrom, masked)
		summary_mat = rbind(summary_mat, dat_temp)
	}
}

#
write.csv(summary_mat, file = "percentages_masked", header=FALSE)
		
