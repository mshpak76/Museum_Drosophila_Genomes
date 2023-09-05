setwd("/raid10/mshpak/MuseumFlies/Counts")

#calculate frequency of inversion-associated alleles at loci linked to common Drosophila inversions

Counts_FR_2L = read.csv("AllCounts_Lyon_Chr2L.txt", header=FALSE, sep='\t')
Counts_Lund_2L = read.csv("AllCounts_ModernLund_Chr2L.txt", header=FALSE, sep='\t')
Counts_1800_2L = read.csv("AllCounts_Nineteenth_Chr2L.txt", header=FALSE, sep='\t')
Counts_1900_2L = read.csv("AllCounts_Twentieth_Chr2L.txt", header=FALSE, sep='\t')

BreakPoints = c(2225744, 13154180)

sum(Counts_FR_2L[BreakPoints[1],])
#182
sum(Counts_FR_2L[BreakPoints[2],])
#50

sum(Counts_Lund_2L[BreakPoints[1],])
#143
sum(Counts_Lund_2L[BreakPoints[2],])
#98

sum(Counts_1800_2L[BreakPoints[1],])
#13
sum(Counts_1800_2L[BreakPoints[2],])
#7

sum(Counts_1900_2L[BreakPoints[1],])
#16
sum(Counts_1900_2L[BreakPoints[2],])
#5


setwd("/raid10/mshpak/MuseumFlies/Inversions/Kapun")

Inv_2Lt = read.csv("In2Lt_fixeddiff_kapun.txt", header=FALSE, sep='\t')

Kapun_Coords_2Lt = Inv_2Lt[,3]
New_Coords_2Lt = Kapun_Coords_2Lt - 1

setwd("/raid10/mshpak/MuseumFlies/Counts")

FR_2L = Counts_FR_2L[Kapun_Coords_2Lt,]
x = rowSums(FR_2L)
# max, min, mean, median = 192, 154, 184.5, 187

Lund20_2L = Counts_Lund_2L[Kapun_Coords_2Lt,]
# max, min, mean, median: 175, 119, 152, 155

Lund_18_2L = Counts_1800_2L[Kapun_Coords_2Lt,]
# 13, 8, 11.3, 12

Lund_19_2L = Counts_1900_2L[Kapun_Coords_2Lt,]
# max, min, mean, median: 16, 12, 14.9, 15.5


Counts_Lund_2R = read.csv("AllCounts_ModernLund_Chr2R.txt", header=FALSE, sep='\t')
Counts_1800_2R = read.csv("AllCounts_Nineteenth_Chr2R.txt", header=FALSE, sep='\t')
Counts_1900_2R = read.csv("AllCounts_Twentieth_Chr2R.txt", header=FALSE, sep='\t')

setwd("/raid10/mshpak/MuseumFlies/Inversions/Kapun")

Inv_2RNs = read.csv("In2RNs_fixeddiff_kapun.txt", header=FALSE, sep='\t')
Lund20_2R = Counts_Lund_2L[Kapun_Coords_2RNs,]
Kapun_Coords_2RNs = Inv_2RNs[,3]

Lund20_2R = Counts_Lund_2R[Kapun_Coords_2RNs,]
#max, min, mean, median 170.50, 0, 139.81, 145.86

Lund19_2R = Counts_1900_2R[Kapun_Coords_2RNs,]
#max, min, mean, median 20, 4, 12.9, 10

Lund18_2R = Counts_1800_2R[Kapun_Coords_2RNs,]
#max, min, mean, median 14, 6, 11.58, 12

setwd("/raid10/mshpak/MuseumFlies/Counts")

Counts_Lund_3R = read.csv("AllCounts_ModernLund_Chr3R.txt", header=FALSE, sep='\t')
Counts_1800_3R = read.csv("AllCounts_Nineteenth_Chr3R.txt", header=FALSE, sep='\t')
Counts_1900_3R = read.csv("AllCounts_Twentieth_Chr3R.txt", header=FALSE, sep='\t')

Counts_Lund_3L = read.csv("AllCounts_ModernLund_Chr3L.txt", header=FALSE, sep='\t')
Counts_1800_3L = read.csv("AllCounts_Nineteenth_Chr3L.txt", header=FALSE, sep='\t')
Counts_1900_3L = read.csv("AllCounts_Twentieth_Chr3L.txt", header=FALSE, sep='\t')


setwd("/raid10/mshpak/MuseumFlies/Inversions/Kapun")
Inv_3LP = read.csv("In3LP_fixeddiff_kapun.txt", header=FALSE, sep='\t')
Kapun_Coords_3LP = Inv_3LP[,3]

Lund20_3L = Counts_Lund_3L[Kapun_Coords_3LP,]
# max, min, mean, median 166.31, 0, 140.85, 147.86

Lund19_3L = Counts_1900_3L[Kapun_Coords_3LP,]
# max, min, mean, median 13, 6, 10.51, 10

Lund18_3L = Counts_1800_3L[Kapun_Coords_3LP,]
# max, min, mean, median 11, 6, 10.29, 11

Inv_3RC = read.csv("In3RC_fixeddiff_kapun.txt", header=FALSE, sep='\t')
Inv_3RK = read.csv("In3RK_fixeddiff_kapun.txt", header=FALSE, sep='\t')
Inv_3RMo = read.csv("In3RMo_fixeddiff_kapun.txt", header=FALSE, sep='\t')
Inv_3RPayne = read.csv("In3RPayne_fixeddiff_kapun.txt", header=FALSE, sep='\t')


Kapun_Coords_3RC = Inv_3RC[,3]

Lund20_3RC = Counts_Lund_3R[Kapun_Coords_3RC,]
#max, min, mean, median 170.67, 0, 145.82, 149.64

Lund19_3RC = Counts_1900_3R[Kapun_Coords_3RC,]
#max, min, mean, median 23, 9, 14.19, 14

Lund18_3RC = Counts_1800_3R[Kapun_Coords_3RC,]
#max, min, mena, median 14, 6, 11.33, 11

Kapun_Coords_3RMo = Inv_3RMo[,3]

Lund20_3RM = Counts_Lund_3R[Kapun_Coords_3RMo,]
#max, min, mean, median 169.69, 0, 144.76, 148.41

Lund19_3RM = Counts_1900_3R[Kapun_Coords_3RMo,]
# max, min, mean, median 23, 7, 13.6, 12

Lund18_3RM = Counts_1800_3R[Kapun_Coords_3RMo,]
#max, min, mena, median 14, 8, 11.71, 12

Kapun_Coords_3RPayne = Inv_3RPayne[,3]

Lund20_3RP = Counts_Lund_3R[Kapun_Coords_3RPayne,]
# max, min, mean, median 163.04, 0, 144.06, 150.72

Lund19_3RP = Counts_1900_3R[Kapun_Coords_3RPayne,]
# max, min, mean, median 18, 12, 15.842, 15

Lund18_3RP = Counts_1800_3R[Kapun_Coords_3RPayne,]
# max, min, mean, meidna 13, 6, 10.68, 10






