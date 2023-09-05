# this version computes pairwise Hamming distances based on sites that are not NA in both genotypes from each pair, as opposed to averaging
# over all polymorphic sites

AllDat=read.csv("All_Poly_Genotype_Matrix_ChrX_RemoveTelomere_For_PCA.csv",header=TRUE)

old_names = names(AllDat)

m18_names = c('M18_4','M18_5','M18_6','M18_10', 'M18_11','M18_12','M18_13','M18_9','M18_25')
m19_names = c('M1','M3','M7','M8','M14','M15','M16','M17','M21','M18','M19','M20','M22','M23','M24')
eg_names = old_names[25:34]
neth_names = old_names[50:59]
#stockholm:69-78
stock_names = c('St1','St2','St3','St4','St5','St6','St7','St8','St9','St10')
#italy: 95-104
italy_names = c('It1','It2','It3','It4','It5','It6','It7','It8','It9','It10')
fr_names = old_names[110:119]

NewDat = cbind(AllDat[,1:24], cbind(AllDat[,25:34], cbind(AllDat[,50:59], cbind(AllDat[,69:78], cbind(AllDat[,95:104],AllDat[,110:119])))))
colnames(NewDat) = c(m18_names, m19_names, eg_names, neth_names, stock_names, italy_names, fr_names)

Hamming = function(x,y){
	get_dist = sum(x != y, na.rm=TRUE)
	return(get_dist)
}
	
	
#as a fraction of non-na sites
Hamming_Proportion = function(x,y){
	ham_counts = Hamming(x,y)
	non_na_sites = length(which(is.na(x)==FALSE & is.na(y)==FALSE))
	return(ham_counts/non_na_sites)
}


# use this - count all sites to determine Hamming distance fraction (otherwise distance is inflated)
#Hamming_Proportion = function(x,y){
#	ham_counts = Hamming(x,y)
#	#non_na_sites = length(which(is.na(x)==FALSE & is.na(y)==FALSE))
#	return(ham_counts/length(x))
#}

nsamples = dim(NewDat)[2]
Dist_Mat = matrix(0, nsamples,nsamples)
for(i in 1:nsamples){
	Dist_Mat[i,i] = 0
	for(j in 1:i){
		Dist_Mat[i,j] = Hamming_Proportion(NewDat[,i],NewDat[,j])
		Dist_Mat[j,i] = Dist_Mat[i,j]
	}
}

colnames(Dist_Mat) = colnames(NewDat)
rownames(Dist_Mat) = colnames(NewDat)

write.csv(Dist_Mat, file="Distance_Matrices_AllPolyLoci_NoNA.csv",row.names=TRUE)

#remove NA
CleanDat = na.omit(NewDat)
nsamples = dim(CleanDat)[2]
Clean_Dist_Mat = matrix(0, nsamples, nsamples)
for(i in 1:nsamples){
	Clean_Dist_Mat[i,i]=0
	for(j in 1:i){
		Clean_Dist_Mat[i,j] = Hamming_Proportion(CleanDat[,i], CleanDat[,j])
		Clean_Dist_Mat[j,i] = Clean_Dist_Mat[i,j]
	}
}


colnames(Clean_Dist_Mat) = colnames(NewDat)
rownames(Clean_Dist_Mat) = colnames(NewDat)

write.csv(Clean_Dist_Mat, file="Distance_Matrices_NoNa.csv",row.names=TRUE)

New_Mat = read.csv("Distance_Matrices_NoNa.csv", header=TRUE, sep=',')
New_Mat = New_Mat[,2:75]
names_1933 = c("M19_1","M19_3","M19_7","M19_8","M19_14","M19_15","M19_16","M19_17","M19_21","M19_18","M19_19","M19_20","M19_22","M19_23","M19_24")
names(New_Mat)[10:24]=names_1933


x = hclust(as.dist(New_Mat))
pdf("Dist_ChromX.pdf")
plot(x, hang = -1, xlab="", sub="")
dev.off()

pdf("Dist_ChromX.pdf")
plot(hclust(as.dist(New_Mat)))
dev.off()


# museum fly comparison: 1800's, 1900's, all modern
#1800
Diffs_1800 = matrix(0,9,3)
for(i in 1:9){
	vtemp1800 = Clean_Dist_Mat[i,1:9][-i]
	mean_temp_1800 = mean(vtemp1800)
	Diffs_1800[i,1]=mean_temp_1800
	vtemp1933 = Clean_Dist_Mat[i,10:24]
	mean_temp_1933 = mean(vtemp1933)
	Diffs_1800[i,2]=mean_temp_1933
	vtemp2000 = Clean_Dist_Mat[i,25:74]
	mean_temp_2000 = mean(vtemp2000)
	Diffs_1800[i,3] = mean_temp_2000
}

#1933
Diffs_1933 = matrix(0,15,3)
for(i in 1:15){
	vtemp1800 = Clean_Dist_Mat[i+9,1:9]
	mean_temp_1800 = mean(vtemp1800)
	Diffs_1933[i,1] = mean_temp_1800
	#print("1800 done")
	vtemp1933 = Clean_Dist_Mat[i+9,10:24][-(i+9)]
	mean_temp_1933 = mean(vtemp1933)
	Diffs_1933[i,2] = mean_temp_1933
	#print("1933 done")1
	vtemp2000 = Clean_Dist_Mat[i+9,25:74]
	mean_temp_2000 = mean(vtemp2000)
	Diffs_1933[i,3] = mean_temp_2000
}

Diffs_Modern = c()
for(i in 1:50){
	vtemp2000 = Clean_Dist_Mat[i+24,25:74][-(i+24)]
	Diffs_Modern = c(Diffs_Modern, mean(vtemp2000))
}

#mean Diffs_Modern (among 21st century samples) : 0.14066
#mean 1800s/1800s: 0.1013
#mean 1800s/1933: 0.1221
#mean 1800s/modern: 0.1295
#mean 1933/1933: 0.1060
#mean 1933/modern: 0.13947
#mean modern/modern: 0.14066
	
	
