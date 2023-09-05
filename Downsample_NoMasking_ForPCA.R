#Perform principal components analysis on samples
#Unrelated individuals only, no masking 

# X chromosome only, remove telomere regions with low recombination
# only select Lund flies from 1800s, downsample the rest to match sample size

setwd("/raid10/mshpak/MuseumFlies/Counts/Unmasked_Counts")

Mus1800 = read.csv("AllCounts_Nineteenth_LundOnly_ChrX.txt",header=FALSE,sep='\t')
s1800 = rowSums(Mus1800)

Mus1933 = read.csv("AllCounts_Twentieth_ChrX.txt",header=FALSE,sep='\t')
s1933 = rowSums(Mus1933)

#Mus_2L = read.csv("AllCounts_Museum_Chr2L.txt",header=FALSE,sep='\t')
#Mus_2R = read.csv("AllCounts_Museum_Chr2R.txt",header=FALSE,sep='\t')
#Mus_3L = read.csv("AllCounts_Museum_Chr3L.txt",header=FALSE,sep='\t')
#Mus_3R = read.csv("AllCounts_Museum_Chr3R.txt",header=FALSE,sep='\t')

Mus = Mus1800+Mus1933

setwd("/raid10/mshpak/MuseumFlies/Counts")

FR = read.csv("AllCounts_Lyon_ChrX.txt",header=FALSE,sep='\t')
sfr = rowSums(FR)


Eg = read.csv("AllCounts_Egypt_ChrX.txt",header=FALSE,sep='\t')
seg = rowSums(Eg)

Neth = read.csv("AllCounts_Netherlands_NoInv_ChrX.txt",header=FALSE,sep='\t')
sneth = rowSums(Neth)

Stock= read.csv("AllCounts_Stockholm_ChrX.txt",header=FALSE,sep='\t')
sstock = rowSums(Stock)

Ital = read.csv("AllCounts_Italy_NoInv_ChrX.txt",header=FALSE,sep='\t')
sital = rowSums(Ital)

to_keep = which(s1800 > 3 & s1800 > 3 & sfr > 3 & seg > 3 & sneth > 3 & sstock > 3 & sital > 3)
to_keep = to_keep[which(to_keep >= 2300000 & to_keep <= 21400000)]

#
All_Dat = Mus + FR + Eg + Neth + Stock + Ital

#rm(list = c('Ital', 'Stock', 'Eg', 'FR', 'Neth', 'Mus', 'Mus1800','Mus1933'))

find_max = function(v){
	tmax = max(v)
	return(which(v==tmax)[1])
}

find_second = function(v){
	v[which(v==max(v))][1]=-1
	tmax=max(v)
	return(which(v==tmax)[1])
}


All_Sum = rowSums(All_Dat)
max_pos = apply(All_Dat, 1, find_max)
All_Freq = All_Dat[cbind(seq_len(nrow(All_Dat)),max_pos)]/All_Sum
#require at least 2 minor alleles
all_poly = which(All_Freq > 0.0357 & All_Freq < 0.964 & All_Sum > 27)
all_poly = all_poly[which(all_poly > 2300000 & all_poly < 21400000)]
all_poly = intersect(all_poly, to_keep)

All_Poly_Mat = All_Dat[all_poly,]

#identify major and minor alleles
#positions of most frequent allele in some population
max_pos_poly = apply(All_Poly_Mat, 1, find_max)
minor_pos_poly = apply(All_Poly_Mat, 1, find_second)

#create haploid genotypes from count file for each individual
#haploid: all EG, N, H4,H5,H6, H12 (1800's), H3, H7 (1933)
#diploid: FR, IT, Stockholm, H9, H10, H11, H13, H25 (1800's), H1, H8, H14-H24 (1933)

rm(list = c('All_Dat', 'All_Poly_Mat', 'Ital', 'Stock', 'Eg', 'FR', 'Neth', 'Mus'))



#HAPLOID
#setwd("/raid10/mshpak/MuseumFlies/IndivCounts/Chrom_X/Haploid")

#Hap1800 = c('H4','H5','H6','H10','H11','H12','H13')
# Lund Only
#Hap1800 = c('H4','H6','H10','H11','H12','H13')
Hap1800 = c('H5','H11','H13')
#Hap1933 = c('H1','H3','H7','H8','H14','H15','H16','H17','H21')
Hap1933 = c('H8','H15','H16')
Netherlands = c('N01','N02','N03','N04','N07','N10','N11','N13','N14','N15','N16','N17','N18','N19','N22','N23','N25','N29','N30')
Netherlands = Netherlands[1:5]
Egypt = c('EG12N','EG15N','EG16N','EG25N','EG26N','EG28N','EG29N','EG31N','EG33N','EG34N','EG35N','EG36N','EG38N','EG44N','EG46N','EG48N','EG49N','EG50N','EG53N','EG54N','EG55N','EG57N','EG58N','EG59N','EG65N','EG70N','EG74N','EG75N','EG76N')
#downsample
Egypt = Egypt[1:5]

#Chroms = c("Chr2L","Chr2R","Chr3L","Chr3R")
Chroms = c("ChrX")

make_genotype = function(maj_allele, min_allele){
	if(is.na(maj_allele) | is.na(min_allele)){
		return(NA)
	}
	else if(maj_allele > 0){
		return(1)
	}else if(min_allele > 0){
		return(0)
	}else{
		return(NA)
	}
}

setwd("/raid10/mshpak/MuseumFlies/IndivCounts/MuseumIndivs")

# also need counts: add rather than take cbind

Mat_Hap_1800 = c() #create a matrix where each column is a haploid genotype
Total_Hap_1800 = matrix(0,22422827,4)[all_poly,]
for(hapname_1800 in Hap1800){
	Sample_Mat = c()
	for(chrom in Chroms){
	        New_Mat = read.csv(paste(c("IndivCounts_",hapname_1800,"_NoInv_",chrom,".txt"),collapse=""), header=FALSE, sep='\t')
	        Sample_Mat = rbind(Sample_Mat, New_Mat)
	}
	Mat_Keep = Sample_Mat[all_poly,]
	Total_Hap_1800 = Total_Hap_1800 + Mat_Keep
	major_alleles = Mat_Keep[cbind(seq_len(nrow(Mat_Keep)), max_pos_poly)]
	minor_alleles = Mat_Keep[cbind(seq_len(nrow(Mat_Keep)), minor_pos_poly)]
	gtypes = mapply(make_genotype, major_alleles, minor_alleles)
	Mat_Hap_1800 = cbind(Mat_Hap_1800, gtypes)
}

colnames(Mat_Hap_1800) = Hap1800
#SMat_1800 = rowSums(Total_1800)



Mat_Hap_1933 = c() #create a matrix where each column is a haploid genotype
Total_Hap_1933 = matrix(0,22422827,4)[all_poly,]
for(hapname_1933 in Hap1933){
	Sample_Mat = c()
	for(chrom in Chroms){
	        New_Mat = read.csv(paste(c("IndivCounts_",hapname_1933,"_NoInv_",chrom,".txt"),collapse=""), header=FALSE, sep='\t')
	        Sample_Mat = rbind(Sample_Mat, New_Mat)
	}
	Mat_Keep = Sample_Mat[all_poly,]
	Total_Hap_1933 = Total_Hap_1933 + Mat_Keep
	major_alleles = Mat_Keep[cbind(seq_len(nrow(Mat_Keep)), max_pos_poly)]
	minor_alleles = Mat_Keep[cbind(seq_len(nrow(Mat_Keep)), minor_pos_poly)]
	gtypes = mapply(make_genotype, major_alleles, minor_alleles)
	Mat_Hap_1933 = cbind(Mat_Hap_1933, gtypes)
}
colnames(Mat_Hap_1933) = Hap1933

setwd("/raid10/mshpak/MuseumFlies/IndivCounts/Chrom_X/Haploid")

Neth_Hap = c() #create a matrix where each column is a haploid genotype
Total_Neth = matrix(0,22422827,4)[all_poly,]
for(hapname_neth in Netherlands){
	Sample_Mat = c()
	for(chrom in Chroms){
	        New_Mat = read.csv(paste(c("IndivCounts_",hapname_neth,"_NoInv_",chrom,".txt"),collapse=""), header=FALSE, sep='\t')
	        Sample_Mat = rbind(Sample_Mat, New_Mat)
	}
	Mat_Keep = Sample_Mat[all_poly,]
	Total_Neth = Total_Neth + Mat_Keep
	major_alleles = Mat_Keep[cbind(seq_len(nrow(Mat_Keep)), max_pos_poly)]
	minor_alleles = Mat_Keep[cbind(seq_len(nrow(Mat_Keep)), minor_pos_poly)]
	gtypes = mapply(make_genotype, major_alleles, minor_alleles)
	Neth_Hap = cbind(Neth_Hap, gtypes)
}
colnames(Neth_Hap) = Netherlands
SNeth = rowSums(Total_Neth)

Egypt_Hap = c() #create a matrix where each column is a haploid genotype
Total_Egypt = matrix(0,22422827,4)[all_poly,]
for(hapname_egypt in Egypt){
	Sample_Mat = c()
	for(chrom in Chroms){
	        New_Mat = read.csv(paste(c("IndivCounts_",hapname_egypt,"_NoInv_",chrom,".txt"),collapse=""), header=FALSE, sep='\t')
	        Sample_Mat = rbind(Sample_Mat, New_Mat)
	}
	Mat_Keep = Sample_Mat[all_poly,]
	Total_Egypt = Total_Egypt + Mat_Keep
	major_alleles = Mat_Keep[cbind(seq_len(nrow(Mat_Keep)), max_pos_poly)]
	minor_alleles = Mat_Keep[cbind(seq_len(nrow(Mat_Keep)), minor_pos_poly)]
	gtypes = mapply(make_genotype, major_alleles, minor_alleles)
	Egypt_Hap = cbind(Egypt_Hap, gtypes)
}
colnames(Egypt_Hap) = Egypt
SEgypt = rowSums(Total_Egypt)
#
	    
	        
	


#DIPLOID
#setwd("/raid10/mshpak/MuseumFlies/IndivCounts/Chrom_X/Diploid")

Dip1800 = c('H9','H25')
#Dip1933 = c('H18','H19','H20','H22','H23','H24')
Dip1933 = c('H20','H23')

Italy = c('SRR5762776','SRR5762780','SRR5762792','SRR5762793','SRR5762794','SRR5762795','SRR5762796','SRR5762797','SRR5762798','SRR5762800','SRR5762801','SRR5762802','SRR5762803','SRR5762804','SRR5762809','SRR5762810')
Italy = Italy[1:5]

Stockholm =c('SRR5762771','SRR5762772','SRR5762773','SRR5762774','SRR5762775','SRR5762777','SRR5762778','SRR5762779','SRR5762781','SRR5762782','SRR5762783',
'SRR5762784','SRR5762785','SRR5762786','SRR5762786','SRR5762787','SRR5762788','SRR5762789','SRR5762790','SRR5762791','SRR5762799','SRR5762805',
'SRR5762806','SRR5762807','SRR5762808','SRR5762811','SRR5762812','SRR6765735','SRR6765736')
#Stockholm = Stockholm[1:25]
Stockholm = Stockholm[1:5]

#remove FR11N, 70N
France = c('FR106N','FR109N','FR112N','FR115N','FR11N','FR126N','FR12N','FR147N','FR14','FR151','FR152N',
'FR153N','FR157N','FR158N','FR162N','FR164N','FR169N','FR16N','FR173N','FR180','FR186N','FR188N','FR19N','FR200N',
'FR207','FR208N','FR213N','FR216N','FR217','FR219N','FR222N','FR225N','FR229','FR230N','FR231N','FR232N',
'FR235N','FR236N','FR238N','FR240N','FR248N','FR252N','FR257N','FR260N','FR261N','FR263N','FR264N','FR269N','FR26N','FR276N',
'FR279N','FR284N','FR288N','FR28N','FR293N','FR296N','FR299N','FR2N','FR302N','FR304N','FR305N','FR310','FR312Nsites',
'FR313N','FR319N','FR320N','FR323N','FR326N','FR32N','FR340N','FR345N','FR348N','FR34N','FR357N','FR360N','FR361',
'FR364N','FR370N','FR37N','FR48N','FR54N','FR59N','FR5N','FR60N','FR73N','FR89N','FR91N','FR94N')
#France = France[1:25]
France = France[1:5]

#sample haploid genotype from diploid data (for consistency with haploid genotyping of other genomes)
# this also works for "mixed" haploid/diploid genotypes such as France, e.g. 1 for IBD regions, 0/1/2 for non IBD
make_diploid_genotype = function(maj_allele, min_allele){
	if(is.na(maj_allele) | is.na(min_allele)){
		return(NA)
	}
	else if(maj_allele > 0 & min_allele == 0){
		return(1)
	}
	else if(maj_allele > 0 & min_allele > 0){
		return(sample(0:1,1))
	}
	else if(maj_allele == 0 & min_allele > 0){
		return(0)
	}
	else{
		return(NA)
	}
}

#alternative: if not maj_allele, set to 0 rather than to NA


Chroms = c("ChrX")

setwd("/raid10/mshpak/MuseumFlies/IndivCounts/MuseumIndivs")

Mat_Dip_1800 = c() #create a matrix where each column is a haploid genotype
Total_Dip_1800 = matrix(0,22422827,4)[all_poly,]
for(dipname_1800 in Dip1800){
	Sample_Mat = c()
	for(chrom in Chroms){
	        New_Mat = read.csv(paste(c("IndivCounts_",dipname_1800,"_NoInv_",chrom,".txt"),collapse=""), header=FALSE, sep='\t')
	        Sample_Mat = rbind(Sample_Mat, New_Mat)
	}
	Mat_Keep = Sample_Mat[all_poly,]
	Total_Dip_1800 = Total_Dip_1800 + Mat_Keep
	major_alleles = Mat_Keep[cbind(seq_len(nrow(Mat_Keep)), max_pos_poly)]
	minor_alleles = Mat_Keep[cbind(seq_len(nrow(Mat_Keep)), minor_pos_poly)]
	gtypes = mapply(make_diploid_genotype, major_alleles, minor_alleles)
	Mat_Dip_1800 = cbind(Mat_Dip_1800, gtypes)
}
colnames(Mat_Dip_1800) = Dip1800
Mat_1800 = cbind(Mat_Hap_1800, Mat_Dip_1800)
Total_1800 = Total_Hap_1800 + Total_Dip_1800
SMat_1800 = rowSums(Total_1800)


Mat_Dip_1933 = c() #create a matrix where each column is a haploid genotype
Total_Dip_1933 = matrix(0,22422827,4)[all_poly,]
for(dipname_1933 in Dip1933){
	Sample_Mat = c()
	for(chrom in Chroms){
	        New_Mat = read.csv(paste(c("IndivCounts_",dipname_1933,"_NoInv_",chrom,".txt"),collapse=""), header=FALSE, sep='\t')
	        Sample_Mat = rbind(Sample_Mat, New_Mat)
	}
	Mat_Keep = Sample_Mat[all_poly,]
	Total_Dip_1933 = Total_Dip_1933 + Mat_Keep
	major_alleles = Mat_Keep[cbind(seq_len(nrow(Mat_Keep)), max_pos_poly)]
	minor_alleles = Mat_Keep[cbind(seq_len(nrow(Mat_Keep)), minor_pos_poly)]
	gtypes = mapply(make_diploid_genotype, major_alleles, minor_alleles)
	Mat_Dip_1933 = cbind(Mat_Dip_1933, gtypes)
}
colnames(Mat_Dip_1933) = Dip1933

Mat_1933 = cbind(Mat_Hap_1933, Mat_Dip_1933)
Total_1933 = Total_Hap_1933 + Total_Dip_1933
SMat_1933 = rowSums(Total_1933)

setwd("/raid10/mshpak/MuseumFlies/IndivCounts/Chrom_X/Diploid")

#sample single haploid genotype for consistency with other data sets
Italy_Dip = c() #create a matrix where each column is a haploid genotype
Total_Italy = matrix(0,22422827,4)[all_poly,]
for(dipname_italy in Italy){
	Sample_Mat = c()
	for(chrom in Chroms){
	        New_Mat = read.csv(paste(c("IndivCounts_",dipname_italy,"_NoInv_",chrom,".txt"),collapse=""), header=FALSE, sep='\t')
	        Sample_Mat = rbind(Sample_Mat, New_Mat)
	}
	Mat_Keep = Sample_Mat[all_poly,]
	Total_Italy = Total_Italy + Mat_Keep
	major_alleles = Mat_Keep[cbind(seq_len(nrow(Mat_Keep)), max_pos_poly)]
	minor_alleles = Mat_Keep[cbind(seq_len(nrow(Mat_Keep)), minor_pos_poly)]
	gtypes = mapply(make_diploid_genotype, major_alleles, minor_alleles)
	Italy_Dip = cbind(Italy_Dip, gtypes)
}
colnames(Italy_Dip) = Italy
SItaly = rowSums(Total_Italy)


Stockholm_Dip = c() #create a matrix where each column is a haploid genotype
Total_Stock = matrix(0,22422827,4)[all_poly,]
for(dipname_stock in Stockholm){
	Sample_Mat = c()
	for(chrom in Chroms){
	        New_Mat = read.csv(paste(c("IndivCounts_",dipname_stock,"_NoInv_",chrom,".txt"),collapse=""), header=FALSE, sep='\t')
	        Sample_Mat = rbind(Sample_Mat, New_Mat)
	}
	Mat_Keep = Sample_Mat[all_poly,]
	Total_Stock = Total_Stock + Mat_Keep
	major_alleles = Mat_Keep[cbind(seq_len(nrow(Mat_Keep)), max_pos_poly)]
	minor_alleles = Mat_Keep[cbind(seq_len(nrow(Mat_Keep)), minor_pos_poly)]
	gtypes = mapply(make_diploid_genotype, major_alleles, minor_alleles)
	Stockholm_Dip = cbind(Stockholm_Dip, gtypes)
}
colnames(Stockholm_Dip) = Stockholm

SStock = rowSums(Total_Stock)


France_Dip = c() #create a matrix where each column is a haploid genotype
Total_France = matrix(0,22422827,4)[all_poly,]
for(dipname_france in France){
	Sample_Mat = c()
	for(chrom in Chroms){
	        New_Mat = read.csv(paste(c("IndivCounts_",dipname_france,"_NoInv_",chrom,".txt"),collapse=""), header=FALSE, sep='\t')
	        Sample_Mat = rbind(Sample_Mat, New_Mat)
	}
	Mat_Keep = Sample_Mat[all_poly,]
	Total_France = Total_France + Mat_Keep
	major_alleles = Mat_Keep[cbind(seq_len(nrow(Mat_Keep)), max_pos_poly)]
	minor_alleles = Mat_Keep[cbind(seq_len(nrow(Mat_Keep)), minor_pos_poly)]
	gtypes = mapply(make_diploid_genotype, major_alleles, minor_alleles)
	France_Dip = cbind(France_Dip, gtypes)
}
colnames(France_Dip) = France

SFrance = rowSums(Total_France)

#merge files for PCA
#
Full_Mat = cbind(Mat_1800, Mat_1933)
Full_Mat = cbind(Full_Mat, cbind(Egypt_Hap,Neth_Hap))
Full_Mat = cbind(Full_Mat, cbind(Stockholm_Dip, Italy_Dip))
Full_Mat = cbind(Full_Mat, France_Dip)

# make sure downsampling sites have 4 or more alleles
to_sample = which(SFrance > 3 & SEgypt > 3 & SStock > 3 & SMat_1933 > 3 & SMat_1800 > 3 & SNeth > 3)
Full_Mat = Full_Mat[to_sample,]

#rm(list=c('France_Dip','Egypt_Hap','Stockholm_Dip','Italy_Dip','Mat_Dip_1800','Mat_Hap_1800','Mat_Hap_1933','Mat_Dip_1933'))

#rename_cols = c('H4','H5','H6','H10','H11','H12','H13','H9','H25','H1','H3','H7','H8','H14','H15','H16','H17','H21')
rename_cols = c('H5','H11','H13','H9','H25','H8','H15','H16','H20','H23')
rename_cols = c(rename_cols, Egypt, Netherlands, Stockholm, Italy, France) 

#temp_names = colnames(Full_Mat)
#temp_names = c(rename_cols, temp_names[19:134])

colnames(Full_Mat) = rename_cols

setwd("/raid10/mshpak/MuseumFlies/Counts/Downsample/Downsample_To_Four")
write.table(Full_Mat, file = "NoMask_Poly_Genotype_Matrix_ChrX_RemoveTelomere_For_PCA.csv", row.names=FALSE, col.names=TRUE, sep='\t')

#
Freqmat = read.csv("NoMask_Poly_Genotype_Matrix_ChrX_RemoveTelomere_For_PCA.csv", header=TRUE, sep='\t')

#determine number of na entries
check_na = function(v){
	return(length(which(is.na(v))))
}


#consider as an option
sample_nas = apply(t(Freqmat), 1, check_na)
which(sample_nas > 200000)
#Freqmat = Freqmat[,-c(1,2,3,10,32,114,117)]
#remove those with most nas

#need to rename for identification
all_names = colnames(Freqmat)

h_1800 = c('H18_5','H18_11','H18_13','H18_9','H18_25')

h_1933 = c('H19_8','H19_15','H19_16','H19_20', 'H19_23')

eg = c('EG12N','EG15N','EG16N','EG25N','EG26N')

net = c('N01','N02','N03','N04','N07')

stock = c('Sk1','Sk2','Sk3','Sk4','Sk5')

ital = c('It1','It2,','It3','It4','It5')

FR = c('FR106N','FR109N','FR112N','FR115N','FR11N')

new_names = c(h_1800, h_1933, eg, net, stock, ital, FR)	

colnames(Freqmat) = new_names

# NA issue for 1800s, downsample to 4
#more_downsample = c(2,4:6,9:16,19:22,25:28,33:40)
#Freqmat_Down = Freqmat[,more_downsample]

newmat = na.omit(Freqmat)
# try alternative approaches to retain rows/columns with large number of NAs, eg. omit problematic columns
pca_test = prcomp(t(newmat), scores = TRUE, cor = TRUE, na.action=na.omit)
#loadings
pca_loadings = pca_test$rotation

x1800=c(22,0,0,15,15)
x1933=rep(1,5)
xeg=rep(2,5)
xnet=rep(3,5)
xstock=rep(4,5)
xit=rep(5,5)
xlyon=rep(6,5)
xall = c(x1800, x1933, xeg, xnet, xstock, xit, xlyon)

scores.df = data.frame(pca_test$x)
pdf("Fig3B_NoMask_Downsample_PCA.pdf")
plot(x = scores.df$PC1, y = scores.df$PC2, pch=xall, xlab="PC 1", ylab="PC 2")
legend(-35,-20,c("Denmark 1850s","Lund 1809", "Passau 1830s", "Lund 1933", "Egypt", "Netherlands", "Stockholm", "Italy", "Lyon"),col=rep("black",5),pch=c(22,0,15,1,2,3,4,5,6))
dev.off()


# scores of samples on 4 PCA axes
scores.df = data.frame(pca_test$x)
plot(x = scores.df$PC1, y = scores.df$PC2)
text(scores.df$PC1, scores.df$PC2, labels = row.names(scores.df))

#pdf("Downsampled_PCA.pdf")
#plot(x = scores.df$PC1, y = scores.df$PC2)
#text(scores.df$PC1, scores.df$PC2, labels = row.names(scores.df))
#dev.off()

library("scatterplot3d")

mus1800 = rep(0,5) # square
mus1933 = rep(1,5) # circle
eg = rep(2,5) # up triangle
neth = rep(3,5) # +
stock = rep(4,5) # x
ital = rep(5,5) # diamond
fran = rep(6,5) # triangle down
all_symbol = c(mus1800,mus1933,eg,neth,stock,ital,fran)


# add figure legend

pdf("NoMask_Downsampled_PCA.pdf")
plot(x = scores.df$PC1, y = scores.df$PC2, pch=all_symbol)
#text(scores.df$PC1, scores.df$PC2, pch = all_symbol)
dev.off()

pdf("NoMask_Downsampled_PCA_3D.pdf")
scatterplot3d(scores.df$PC1, scores.df$PC2, scores.df$PC3, pch = all_symbol)
dev.off()

#to compare indivmats
whichmax = function(v){
   tempmax = max(v)
   if(tempmax > 0){
   	return(which(v==tempmax)[1])
   }
   else{
   	return(NA)
   	}
   }
