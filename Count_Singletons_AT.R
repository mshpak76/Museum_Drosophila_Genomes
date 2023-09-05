setwd("/raid10/mshpak/MuseumFlies/IndivCounts/MuseumIndivs")

# determine the fraction of singletons that are A/T sites (possible enrichment due to cytosine deamination)

#1800: H4, H5, H6, H9, H10-13, H25
#low quality: H4,H5,H12

#1933: H1, H3, H7, H8, H14-H24
# low qualtiy: H3, H7

Samp_1800 = c("H4","H5","H6","H9","H10","H11","H12","H13","H25")
Samp_1933 = c("H1","H3","H7","H8","H14","H15","H16","H17","H18","H19","H20","H21","H22","H23","H24")

Samp_HQ_1800 = c("H9","H10","H11","H12","H13","H25")
Samp_HQ_1933 =c("H1","H8","H14","H15","H16","H17","H18","H19","H20","H21","H22","H23","H24")

Samp_HQ = c(Samp_HQ_1800, Samp_HQ_1933)
Samp_LQ = c("H4","H5","H6","H12","H3","H7")

fnamebase = "IndivCounts_"
suffix_2L = "_NoInv_Chr2L.txt"
suffix_2R = "_NoInv_Chr2R.txt"
suffix_3L = "_NoInv_Chr3L.txt"
suffix_3R = "_NoInv_Chr3R.txt"
suffix_X = "_NoInv_ChrX.txt"

Temp_2L = read.csv("IndivCounts_H4_NoInv_Chr2L.txt",header=FALSE,sep='\t')
Temp_2R = read.csv("IndivCounts_H4_NoInv_Chr2R.txt",header=FALSE,sep='\t')
Temp_3L = read.csv("IndivCounts_H4_NoInv_Chr3L.txt",header=FALSE,sep='\t')
Temp_3R = read.csv("IndivCounts_H4_NoInv_Chr3R.txt",header=FALSE,sep='\t')
Temp_X = read.csv("IndivCounts_H4_NoInv_ChrX.txt",header=FALSE,sep='\t')
Temp = rbind(Temp_2L, rbind(Temp_2R, rbind(Temp_3L, rbind(Temp_3R, Temp_X))))

rsum = rowSums(Temp)
temp_no_na_sites = which(rsum > 0)

L = dim(Temp)[1]

hq_no_na_sites = seq(1,L,1)
H0 = matrix(0,L,4)
H_HighQual = H0
for(samp in Samp_HQ){
	fnamesamp = paste(fnamebase, samp, sep='')
	fname2L = paste(fnamesamp, suffix_2L, sep='')
	fname2R = paste(fnamesamp, suffix_2R, sep='')
	fname3L= paste(fnamesamp, suffix_3L, sep='')
	fname3R= paste(fnamesamp, suffix_3R, sep='')
	fnameX= paste(fnamesamp, suffix_X, sep='')
	Temp_2L = read.csv(fname2L, header=FALSE,sep='\t')
	Temp_2R = read.csv(fname2R, header=FALSE,sep='\t')
	Temp_3L = read.csv(fname3L, header=FALSE,sep='\t')
	Temp_3R = read.csv(fname3R, header=FALSE,sep='\t')
	Temp_X = read.csv(fnameX, header=FALSE,sep='\t')
	Htemp = rbind(Temp_2L, rbind(Temp_2R, rbind(Temp_3L, rbind(Temp_3R, Temp_X))))
	temp_sum = rowSums(Htemp)
	sites_temp = which(temp_sum > 0)
	hq_no_na_sites = intersect(hq_no_na_sites, sites_temp)
	H_HighQual = H_HighQual + Htemp
}
#

lq_no_na_sites = seq(1,L,1)
H_LowQual = H0
for(samp in Samp_LQ){
	fnamesamp = paste(fnamebase, samp, sep='')
	fname2L = paste(fnamesamp, suffix_2L, sep='')
	fname2R = paste(fnamesamp, suffix_2R, sep='')
	fname3L= paste(fnamesamp, suffix_3L, sep='')
	fname3R= paste(fnamesamp, suffix_3R, sep='')
	fnameX= paste(fnamesamp, suffix_X, sep='')
	Temp_2L = read.csv(fname2L, header=FALSE,sep='\t')
	Temp_2R = read.csv(fname2R, header=FALSE,sep='\t')
	Temp_3L = read.csv(fname3L, header=FALSE,sep='\t')
	Temp_3R = read.csv(fname3R, header=FALSE,sep='\t')
	Temp_X = read.csv(fnameX, header=FALSE,sep='\t')
	Htemp = rbind(Temp_2L, rbind(Temp_2R, rbind(Temp_3L, rbind(Temp_3R, Temp_X))))
	temp_sum = rowSums(Htemp)
	sites_temp = which(temp_sum > 0)
	lq_no_na_sites = intersect(lq_no_na_sites, sites_temp)
	H_LowQual = H_LowQual + Htemp
}

#matrix of non-missing data
# low qual, high qual, all
no_na_sites = intersect(lq_no_na_sites, hq_no_na_sites)

lq_add = rep(NA, length(hq_no_na_sites)-length(lq_no_na_sites))
no_add = rep(NA, length(hq_no_na_sites)-length(no_na_sites))

none = c(no_na_sites, no_add)
lq_only = c(lq_no_na_sites, lq_add)


Non_Na_Mat = rbind(hq_no_na_sites, rbind(lq_only, none))

write.table(H_HighQual, file="AllCount_Genome_ExcludeLowqual.csv", row.names=FALSE, col.names=FALSE, sep='\t')
write.table(H_LowQual, file="AllCount_Genome_Lowqual.csv", row.names=FALSE, col.names=FALSE, sep='\t')
write.table(Non_Na_Mat, file="No_Na_Genome.csv", row.names=FALSE, col.names=FALSE, sep = '\t')
#

# take X chromosome as separate 119029689-22422827:119029689
# allcounts: H_Highqual + H_LowQual

AllDat_HighQual = read.csv("AllCount_Genome_ExcludeLowqual.csv", header=FALSE, sep='\t')
AllDat_LowQual = read.csv("AllCount_Genome_Lowqual.csv", header=FALSE, sep='\t')
No_Missing = read.csv("No_Na_Genome.csv", header=FALSE, sep = '\t')



#to find polymorphic - first identify non-zero sites
poly = function(v){
	return(length(which(v > 0)))
}

H_All = H_HighQual + H_LowQual
rsum_all = rowSums(H_All)
#mean rsum_all 36.7
# those with at least 3 alleles, so that a singleton is a minor allele in some meaningful sense
H_All_Multi = H_All[which(rsum_all > 2),]
# singleton sites
single_all = which(H_All_Multi[,1] == 1 | H_All_Multi[,2] == 1 | H_All_Multi[,3] == 1 | H_All_Multi[,4] == 1)
# singletons as fraction of polymorphic sites
poly_all = apply(H_All_Multi, 1, poly)
psites_all = length(which(poly_all > 1))

sall = intersect(single_all, which(poly_all > 1))
length(sall)/length(which(poly_all > 1))
# 0.4235364

rsum_hq = rowSums(H_HighQual)
# mean rsum_hq = 32.03
# at least 3 alleles to define singleton
H_HQ_Multi = H_HighQual[which(rsum_hq > 2),]
# singleton sites
single_hq = which(H_HQ_Multi[,1] == 1 | H_HQ_Multi[,2] == 1 | H_HQ_Multi[,3] == 1 | H_HQ_Multi[,4] == 1)
# singletons as fraction of polymorphic sites
poly_hq = apply(H_HQ_Multi, 1, poly)
psites_hq = length(which(poly_hq > 1))

shq = intersect(single_hq, which(poly_hq > 1))
length(shq)/length(which(poly_hq > 1))
# 0.3076099
# 36.7/32.03 factor 1.145801, rescaling singleton count to 0.3524495

# low quality sites contribute 27.37% of singletons, despite being 12.73% of alleles 

rsum_lq = rowSums(H_LowQual)
#mean(rsum_lq)
# at least 3 alleles to define singleton
H_LQ_Multi = H_LowQual[which(rsum_lq > 2),]
# singleton sites
single_lq = which(H_LQ_Multi[,1] == 1 | H_LQ_Multi[,2] == 1 | H_LQ_Multi[,3] == 1 | H_LQ_Multi[,4] == 1)
# singletons as fraction of polymorphic sites
poly_lq = apply(H_LQ_Multi, 1, poly)
psites_lq = length(which(poly_lq > 1))

slq = intersect(single_lq, which(poly_lq > 1))
length(slq)/length(which(poly_lq > 1))
#  0.8055777


# check AT content (both overall and singletons)
All_AT = sum(H_All[,1] + H_All[,4])/sum(rowSums(H_All))
# 0.5806435

HQ_AT = sum(H_HighQual[,1] + H_HighQual[,4])/sum(rowSums(H_HighQual))
# 0.5814255

LQ_AT = sum(H_LowQual[,1] + H_LowQual[,4])/sum(rowSums(H_LowQual))
# 0.5752808


# SINGLETONS
H_All_Singletons = H_All_Multi[single_all,]


length(which(H_All_Singletons[,1]==1))+length(which(H_All_Singletons[,4]==1))
# 898062
length(which(H_All_Singletons[,2]==1))+length(which(H_All_Singletons[,3]==1))
# 227886
# AT content = 0.7975273 


H_HQ_Singletons = H_HQ_Multi[single_hq,]

length(which(H_HQ_Singletons[,1]==1))+length(which(H_HQ_Singletons[,4]==1))
# 476601
length(which(H_HQ_Singletons[,2]==1))+length(which(H_HQ_Singletons[,3]==1))
# 167575
# AT content = 0.7398615

# ~ 74-80% of singletons are A/T even though only 58% of overall sites are A/T


##########################################
### NO NA SITES (no missing data)

H_NoNA_All = H_All[no_na_sites,]
H_NoNA_HighQual = H_HighQual[no_na_sites,]
H_NoNA_LowQual = H_LowQual[no_na_sites,]


rsum_all = rowSums(H_NoNA_All)
#mean rsum_all 41.46954
# those with at least 3 alleles, so that a singleton is a minor allele in some meaningful sense
H_All_Multi = H_NoNA_All[which(rsum_all > 2),]
# singleton sites
single_all = which(H_All_Multi[,1] == 1 | H_All_Multi[,2] == 1 | H_All_Multi[,3] == 1 | H_All_Multi[,4] == 1)
# singletons as fraction of polymorphic sites
poly_all = apply(H_All_Multi, 1, poly)
psites_all = length(which(poly_all > 1))
sall = intersect(single_all, which(poly_all > 1))
length(sall)/length(which(poly_all > 1))
# 0.5631823

rsum_hq = rowSums(H_HighQual)
# mean rsum_hq = 32.03115
# at least 3 alleles to define singleton
H_HQ_Multi = H_NoNA_HighQual[which(rsum_hq > 2),]
# singleton sites
single_hq = which(H_HQ_Multi[,1] == 1 | H_HQ_Multi[,2] == 1 | H_HQ_Multi[,3] == 1 | H_HQ_Multi[,4] == 1)
# singletons as fraction of polymorphic sites
poly_hq = apply(H_HQ_Multi, 1, poly)
psites_hq = length(which(poly_hq > 1))

shq = intersect(single_hq, which(poly_hq > 1))
length(shq)/length(which(poly_hq > 1))
# 0.3843134
# mean allele count ratio of all to high qual: 1.294663, correction factor gives 0.4975563
# low quality: 14.47% of sites, 31.76% of singletons

# low quality sites contribute 27.37% of singletons, despite being 12.73% of alleles 

rsum_lq = rowSums(H_NoNA_LowQual)
#mean(rsum_lq)
# at least 3 alleles to define singleton
H_LQ_Multi = H_NoNA_LowQual[which(rsum_lq > 2),]
# singleton sites
single_lq = which(H_LQ_Multi[,1] == 1 | H_LQ_Multi[,2] == 1 | H_LQ_Multi[,3] == 1 | H_LQ_Multi[,4] == 1)
# singletons as fraction of polymorphic sites
poly_lq = apply(H_LQ_Multi, 1, poly)
psites_lq = length(which(poly_lq > 1))

slq = intersect(single_lq, which(poly_lq > 1))
length(slq)/length(which(poly_lq > 1))
# 0.8335462



### AT content for data with no missing sites

All_AT = sum(H_NoNA_All[,1] + H_NoNA_All[,4])/sum(rowSums(H_NoNA_All))
# 0.7151213

HQ_AT = sum(H_NoNA_HighQual[,1] + H_NoNA_HighQual[,4])/sum(rowSums(H_NoNA_HighQual))
# 0.7155148

LQ_AT = sum(H_NoNA_LowQual[,1] + H_NoNA_LowQual[,4])/sum(rowSums(H_NoNA_LowQual))
# 0.7127955


# SINGLETONS
H_NoNA_All_Singletons = H_All_Multi[single_all,]


length(which(H_NoNA_All_Singletons[,1]==1))+length(which(H_NoNA_All_Singletons[,4]==1))
# 143114
length(which(H_NoNA_All_Singletons[,2]==1))+length(which(H_NoNA_All_Singletons[,3]==1))
# 48330
# AT content = 0.7475502 


H_HQ_Singletons = H_HQ_Multi[single_hq,]

length(which(H_HQ_Singletons[,1]==1))+length(which(H_HQ_Singletons[,4]==1))
# 2876
length(which(H_HQ_Singletons[,2]==1))+length(which(H_HQ_Singletons[,3]==1))
# 1028
# AT content = 0.6674844

##### LYON/FR ########
# compare to Lyon flies
setwd("raid10/mshpak/MuseumFlies/Counts")
FR_2L = read.csv("AllCounts_FR_Chr2L.txt", header=FALSE, sep='\t')
FR_2R = read.csv("AllCounts_FR_Chr2R.txt", header=FALSE, sep='\t')
FR_3L = read.csv("AllCounts_FR_Chr3L.txt", header=FALSE, sep='\t')
FR_3R = read.csv("AllCounts_FR_Chr3R.txt", header=FALSE, sep='\t')
FR_X = read.csv("AllCounts_FR_ChrX.txt", header=FALSE, sep='\t')
FR = rbind(FR_2L, rbind(FR_2R, rbind(FR_3L, rbind(FR_3R, FR_X))))



poly = function(v){
	return(length(which(v > 0)))
}

rsum_all = rowSums(FR)
#mean rsum_all 104.368
# those with at least 3 alleles, so that a singleton is a minor allele in some meaningful sense
H_All_Multi = FR[which(rsum_all > 2),]
# singleton sites
single_all = which(H_All_Multi[,1] == 1 | H_All_Multi[,2] == 1 | H_All_Multi[,3] == 1 | H_All_Multi[,4] == 1)
# singletons as fraction of polymorphic sites
poly_all = apply(H_All_Multi, 1, poly)
psites_all = length(which(poly_all > 1))

sall = intersect(single_all, which(poly_all > 1))
length(sall)/length(which(poly_all > 1))
# 0.2515082

# check AT content (both overall and singletons)
All_AT = sum(FR[,1] + FR[,4])/sum(rowSums(FR))
# 0.5712715


# SINGLETONS
H_All_Singletons = H_All_Multi[single_all,]


length(which(H_All_Singletons[,1]==1))+length(which(H_All_Singletons[,4]==1))
# 549063
length(which(H_All_Singletons[,2]==1))+length(which(H_All_Singletons[,3]==1))
# 299180
# ratio of AT among singletons 0.6472945

## NO NA For FR
no_na_fr = which(rsum_all > 0)
no_na = intersect(




##### COMPARE TO LYON No Inversions
setwd("/raid10/mshpak/MuseumFlies/Counts/Inversion_Free/Haploid")

FR_2L = read.csv("AllCounts_FR_NoInv_Chr2L.txt", header=FALSE, sep = '\t')
FR_2R = read.csv("AllCounts_FR_NoInv_Chr2R.txt", header=FALSE, sep = '\t')
FR_3L = read.csv("AllCounts_FR_NoInv_Chr3L.txt", header=FALSE, sep = '\t')
FR_3R = read.csv("AllCounts_FR_NoInv_Chr3R.txt", header=FALSE, sep = '\t')
FR_X = read.csv("AllCounts_FR_NoInv_ChrX.txt", header=FALSE, sep = '\t')

FR = rbind(FR_2L, FR_2R, FR_3L, FR_3R, FR_X)

rsum_all = rowSums(FR)
#mean rsum_all 
# those with at least 3 alleles, so that a singleton is a minor allele in some meaningful sense
H_All_Multi = FR[which(rsum_all > 2),]
# singleton sites
single_all = which(H_All_Multi[,1] == 1 | H_All_Multi[,2] == 1 | H_All_Multi[,3] == 1 | H_All_Multi[,4] == 1)
# singletons as fraction of polymorphic sites
poly_all = apply(H_All_Multi, 1, poly)
psites_all = length(which(poly_all > 1))

sall = intersect(single_all, which(poly_all > 1))
length(sall)/length(which(poly_all > 1))
# 0.2859923

# check AT content (both overall and singletons)
All_AT = sum(FR[,1] + FR[,4])/sum(rowSums(FR))
# 0.5701163


# SINGLETONS
H_All_Singletons = H_All_Multi[single_all,]


length(which(H_All_Singletons[,1]==1))+length(which(H_All_Singletons[,4]==1))
# 477093
length(which(H_All_Singletons[,2]==1))+length(which(H_All_Singletons[,3]==1))
# 263607
# ratio of AT among singletons 0.644111
