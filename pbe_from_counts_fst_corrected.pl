#!/usr/bin/perl -w
#use strict;

#Input - genomic data from exactly three population samples

#Window output file:  
#Three pi and three Dxy values for each window.
#Three FST values calculated from those.
#Window PBS and window PBE values calculated from those.
#SNP Max PBS and SNP Max PBE values for each window.

#SNP output file: (match output of pbs_snp.pl but add PBE column)
#SNP position
#PBE and PBS
#Three FST values (Reynolds code)
# 3 x 4 = 12 site nucleotide counts

my $chr = 'Chr2R';
# also 'ChrX', 'Chr2L', 'Chr3L', 'Chr3R' 'Chr2R'
#my $SNPsPerWindow = 250;
#my $WindowFile = "windows_ZI" . $SNPsPerWindow . "_" . $chr . ".txt";
#my $WindowFile = "windows_10kb_" . $chr . ".txt";
my $WindowFile = "windows_ZI250_" . $chr . ".txt";


#my @pops = ('Nineteenth','ModernLund','Twentieth');
#my @pops = ('Nineteenth','ModernLund','Stockholm');
#my @pops = ('Nineteenth','FR', 'Stockholm');
#my @pops = ('Nineteenth', 'ModernLund', 'FR');
#my @pops = ('Nineteenth','ModernLund', 'FR');
#my @pops = ('Nineteenth', 'ModernLund', 'Lyon');
my @pops = ('Twentieth', 'ModernLund', 'Lyon');
#my @pops = ('Twentieth', 'Stockholm', 'FR');
#my @pops = ('Twentieth','Stockholm','ModernLund');
#my @pops = ('Museum', 'Stockholm', 'FR');
#my @pops = ('Museum', 'Stockholm', 'ModernLund');
#my @pops = ('Museum', 'ModernLund', 'FR'); #combine 1800's and 1933 flies

#Lund samples only for nineteenth
#my @pops = ('NineteenthLund','ModernLund', 'Twentieth');
#my @pops = ('NineteenthLund', 'ModernLund', 'Stockholm');
#my @pops = ('NineteenthLund', 'ModernLund', 'FR'); 

die if (@pops != 3);
my $name = $pops[0] . '_' . $pops[1] . '_' . $pops[2] . 'ZI';
my $OutputFile = $name . '_WinPBE_'. $chr . '.txt';
my $SNPFile = $name . '_SitePBE_'. $chr . '.txt';

my $IncludeInversions = 1;  #0 to exclude arms with inversions, 1 to include

my $SampleMin = 5;  #all pops must have this many alleles called for site PBE/PBS
# change to higher threshold?
# my $SampleMin = 10;
my $MinCov = 100;  #min sites in window meeting the above threshold to report window stats
my $FreqThresh = 0.05;  #avg pop freq needed to report site PBE/PBS in SNPFile (but all sites passing SampleMin are used for window stats)
# possibly make larger to eliminate singletons in reduced museum fly samples - check nucleotide counts per site, i.e. <10 for 1800s
# my $FreqThresh = 0.067;
# i.e. consinder singleton in sample of size 5, frequency 0.2, average freq 0.06666
my $MaxFST = 0.95;  #limit FST to this max value, because we can't use FST=1 in PBE/PBS formula

my $p = 0;
my $q = 0;
my $f = 0;
my $j = 0;
my $s = 0;
my $file = '';
my $file1 = '';
my $FST = 0;
my $kb = -1;
my $window = 0;
my $SpliceNum = 0;
my $SplicePos = 0;
my $SitesAnalyzed = 0;
my $As = 0;
my $Cs = 0;
my $Gs = 0;
my $Ts = 0;
my $AFreqSum = 0;
my $CFreqSum = 0;
my $GFreqSum = 0;
my $TFreqSum = 0;
my $MajorAllele = 0;
my $MinorAllele = 0;
my $WinPBE = 0;
my $WinPBS = 0;
my $SitePBE = 0;
my $SitePBS = 0;
my $MaxPBS = 0;
my $MaxPBE = 0;
my $n = 0;
my $ThreshPassed = 0;
my $AvgFreq = 0;
my $pos = 0;
my $SampleSize1 = 0;
my $SampleSize2 = 0;
my $SampleSize3 = 0;
my $MajorFreq1 = 0;
my $MajorFreq2 = 0;
my $MajorFreq3 = 0;
my $MinorFreq1 = 0;
my $MinorFreq2 = 0;
my $MinorFreq3 = 0;
my $NumA = 0;
my $SharedNum = 0;
my $FracNum = 0;
my $FracDen = 0;
my $frac = 0;
my $WholeNum = 0;
my $WholeDen = 0;
my $DenFracNum = 0;
my $DenFrac = 0;
my $Z = 0;
my $MedianZ = 0;
my $MedianPBS = 0;
my $Num12Sum = 0;
my $Den12Sum = 0;
my $Num13Sum = 0;
my $Den13Sum = 0;
my $Num23Sum = 0;
my $Den23Sum = 0;
my $Dxy12Sum = 0;
my $Dxy13Sum = 0;
my $Dxy23Sum = 0;
my $Pi1Sum = 0;
my $Pi2Sum = 0;
my $Pi3Sum = 0;
my $quantile = 0;
my $site = -1;
my @line = ();
my @Pop1Nucs = ();
my @Pop2Nucs = ();
my @Pop3Nucs = ();
my @AllNucCounts = ();
my @FreqSums = ();
my @WinSiteCounts = ();
my @positions = ();
my @DiffAoA = ();
my @CompAoA = ();
my @SiteCovs = ();
my @PiAoA = ();
my @InputAoA = ();
my @NextAoA = ();
my @DxyAoA = ();
my @FSTAoA = ();
my @pis = ();
my @Dxys = ();
my @FSTs = ();
my @WinPBSs = ();
my @WinPBEs = ();
my @MaxPBSs = ();
my @MaxPBEs = ();
my @NucCounts = ();
my @NucFreqAoA = ();
my @SortedValues = ();
my @WinPBEQuantiles = ();
my @WinPBSQuantiles = ();
my @MaxPBSQuantiles = ();
my @MaxPBEQuantiles = ();
my @WinZs = ();
my @SiteFSTAoA = ();
my @NucCountAoA = ();
my @SiteZs = ();
my @SitePBSs = ();
my @SitePBEs = ();
my @SiteWins = ();


#for ($i = 0; $i < @pops; $i++){
#  push @line, 0;
#}
#for ($i = 0; $i < @pops; $i++){
#  push @DiffAoA, [ @line ];
#  push @CompAoA, [ @line ];
#}

#Get window boundaries
my @WinStarts = ();
my @WinStops = ();
open W, "<$WindowFile" or die "Can't open window file $WindowFile\n";
while (<W>){
  chomp;
  last if m/^$/;
  @line = split;
  push @WinStarts, $line[0];
  push @WinStops, $line[1];
}
close W;

#open count files for each population and get count data site by site
if ($IncludeInversions == 1){
  $file = 'AllCounts_' . $pops[0] . '_' . $chr . '.txt';
  open A, "$file" or die "can not open $file";
  $file = 'AllCounts_' . $pops[1] . '_' . $chr . '.txt';
  open B, "$file" or die "can not open $file";
  $file = 'AllCounts_' . $pops[2] . '_' . $chr . '.txt';
  open C, "$file" or die "can not open $file";
}
else{
  $file = 'AllCounts_' . $pops[0] . '_NoInv_' . $chr . '.txt';
  open A, "$file" or die;
  $file = 'AllCounts_' . $pops[1] . '_NoInv_' . $chr . '.txt';
  open B, "$file" or die;
  $file = 'AllCounts_' . $pops[2] . '_NoInv_' . $chr . '.txt';
  open C, "$file" or die;
}

#open count files for each population and get count data site by site
#if ($IncludeInversions == 1){
#  $file = 'AllCounts_' . $pops[0] . '_' . $chr . '.txt';
#  open A, "<../counts/$file" or die "can not open $file";
#  $file = 'AllCounts_' . $pops[1] . '_' . $chr . '.txt';
#  open B, "<../counts/$file" or die "can not open $file";
#  $file = 'AllCounts_' . $pops[2] . '_' . $chr . '.txt';
#  open C, "<../counts/$file" or die "can not open $file";
#}
#else{
#  $file = 'AllCounts_' . $pops[0] . '_NoInv_' . $chr . '.txt';
#  open A, "<../counts/$file" or die;
#  $file = 'AllCounts_' . $pops[1] . '_NoInv_' . $chr . '.txt';
#  open B, "<../counts/$file" or die;
#  $file = 'AllCounts_' . $pops[2] . '_NoInv_' . $chr . '.txt';
#  open C, "<../counts/$file" or die;
#}

while (<A>){
  chomp;
  @Pop1Nucs = split;
  $_ = (<B>);
  chomp;
  @Pop2Nucs = split;
  $_ = (<C>);
  chomp;
  @Pop3Nucs = split;
  $site++;
  @AllNucCounts = @Pop1Nucs;
  push @AllNucCounts, @Pop2Nucs;
  push @AllNucCounts, @Pop3Nucs;

#if we're at the end of a window, calculate and record window stats
  if ($site > $WinStops[$window]){
    push @WinSiteCounts, $SitesAnalyzed;
    if ($SitesAnalyzed < $MinCov){
      @pis = (-999, -999, -999);
      push @PiAoA, [ @pis ];
      push @DxyAoA, [ @pis ];
    }
    else{
      $Pi1Sum = $Pi1Sum / $SitesAnalyzed;
      $Pi2Sum = $Pi2Sum / $SitesAnalyzed;
      $Pi3Sum = $Pi3Sum / $SitesAnalyzed;
      @pis = ($Pi1Sum, $Pi2Sum, $Pi3Sum);
      push @PiAoA, [ @pis ];
      $Dxy12Sum = $Dxy12Sum / $SitesAnalyzed;
      $Dxy13Sum = $Dxy13Sum / $SitesAnalyzed;
      $Dxy23Sum = $Dxy23Sum / $SitesAnalyzed;
      @Dxys = ($Dxy12Sum, $Dxy13Sum, $Dxy23Sum);
      push @DxyAoA, [ @Dxys ];
    }
    $SitesAnalyzed = 0;
    push @MaxPBSs, $MaxPBS;
    $MaxPBS = 0;
    if (($Den12Sum == 0) || ($Den13Sum == 0) || ($Den23Sum == 0)){
      push @WinPBSs, -999;
      push @WinZs, -999;
      @FSTs = (-999, -999, -999);
      push @FSTAoA, [ @FSTs ];
    }
    else{
      @FSTs = ();
      $FST = $Num12Sum / $Den12Sum;
      if ($FST >= 1){
	$FST = $MaxFST;
      }
      push @FSTs, $FST;
      $FST = $Num13Sum / $Den13Sum;
      if ($FST >= 1){
	$FST = $MaxFST;
      }
      push @FSTs, $FST;
      $FST = $Num23Sum / $Den23Sum;
      if ($FST >= 1){
	$FST = $MaxFST;
      }
      push @FSTs, $FST;
      push @FSTAoA, [ @FSTs ];
      $WinPBS = ((-1 * log (1 - $FSTs[0])) + (-1 * log (1 - $FSTs[1])) - (-1 * log (1 - $FSTs[2]))) / 2;
      push @WinPBSs, $WinPBS;
      $Z = -1 * log (1 - $FSTs[2]);
      push @WinZs, $Z;
    }
    $Pi1Sum = 0;
    $Pi2Sum = 0;
    $Pi3Sum = 0;
    $Dxy12Sum = 0;
    $Dxy13Sum = 0;
    $Dxy23Sum = 0;
    $Num12Sum = 0;
    $Num23Sum = 0;
    $Num13Sum = 0;
    $Den12Sum = 0;
    $Den23Sum = 0;
    $Den13Sum = 0;
    $window++;
    print "got data for $chr, window $window\n";
#    last if ($window == @WinStarts);
  }

#Define major and minor alleles (mask third/fourth alleles)
  $SampleSize1 = $Pop1Nucs[0] + $Pop1Nucs[1] + $Pop1Nucs[2] + $Pop1Nucs[3];
  $SampleSize2 = $Pop2Nucs[0] + $Pop2Nucs[1] + $Pop2Nucs[2] + $Pop2Nucs[3];
  $SampleSize3 = $Pop3Nucs[0] + $Pop3Nucs[1] + $Pop3Nucs[2] + $Pop3Nucs[3];
  next if ($SampleSize1 < $SampleMin);
  next if ($SampleSize2 < $SampleMin);
  next if ($SampleSize3 < $SampleMin);
  $AFreqSum = ($Pop1Nucs[0] / $SampleSize1) + ($Pop2Nucs[0]  / $SampleSize2) + ($Pop3Nucs[0] / $SampleSize3);
  $CFreqSum = ($Pop1Nucs[1] / $SampleSize1) + ($Pop2Nucs[1]  / $SampleSize2) + ($Pop3Nucs[1] / $SampleSize3);
  $GFreqSum = ($Pop1Nucs[2] / $SampleSize1) + ($Pop2Nucs[2]  / $SampleSize2) + ($Pop3Nucs[2] / $SampleSize3);
  $TFreqSum = ($Pop1Nucs[3] / $SampleSize1) + ($Pop2Nucs[3]  / $SampleSize2) + ($Pop3Nucs[3] / $SampleSize3);
  @FreqSums = ($AFreqSum, $CFreqSum, $GFreqSum, $TFreqSum);
  @SortedValues = @FreqSums;
  @SortedValues = sort { $b <=> $a } @SortedValues;
  $MajorAllele = -1;
  for ($j = 0; $j < @FreqSums; $j++){
    if (($FreqSums[$j] == $SortedValues[0]) && ($MajorAllele < -0.5)){
      $MajorAllele = $j;
      next;
    }
    if ($FreqSums[$j] == $SortedValues[1]){
      $MinorAllele = $j;
      next;
    }
  }
  $SampleSize1 = $Pop1Nucs[$MajorAllele] + $Pop1Nucs[$MinorAllele];
  $SampleSize2 = $Pop2Nucs[$MajorAllele] + $Pop2Nucs[$MinorAllele];
  $SampleSize3 = $Pop3Nucs[$MajorAllele] + $Pop3Nucs[$MinorAllele];
  next if ($SampleSize1 < $SampleMin);
  next if ($SampleSize2 < $SampleMin);
  next if ($SampleSize3 < $SampleMin);
  $SitesAnalyzed++;
  
#Pi and Dxy calculations (adding site data to window sums)
  $MajorFreq1 = $Pop1Nucs[$MajorAllele] / $SampleSize1;
  $MajorFreq2 = $Pop2Nucs[$MajorAllele] / $SampleSize2;
  $MajorFreq3 = $Pop3Nucs[$MajorAllele] / $SampleSize3;
  $MinorFreq1 = $Pop1Nucs[$MinorAllele] / $SampleSize1;
  $MinorFreq2 = $Pop2Nucs[$MinorAllele] / $SampleSize2;
  $MinorFreq3 = $Pop3Nucs[$MinorAllele] / $SampleSize3;
  $Pi1Sum += 2 * $MajorFreq1 * $MinorFreq1;
  $Pi2Sum += 2 * $MajorFreq2 * $MinorFreq2;
  $Pi3Sum += 2 * $MajorFreq3 * $MinorFreq3;
  $Dxy12Sum += $MajorFreq1 * $MinorFreq2 + $MajorFreq2 * $MinorFreq1;
  $Dxy13Sum += $MajorFreq1 * $MinorFreq3 + $MajorFreq3 * $MinorFreq1;
  $Dxy23Sum += $MajorFreq3 * $MinorFreq2 + $MajorFreq2 * $MinorFreq3;
  
#Reynolds FST calculations (pops 1 & 2)
  @FSTs = ();
  $SharedNum = $SampleSize1 * ($MinorFreq1 -  ($MinorFreq1 ** 2)) + $SampleSize2 * ($MinorFreq2 - ($MinorFreq2 ** 2));
  # check if NumA needs a factor of 1/2 as per Reynolds et al 1983 formula
  $NumA = ($MinorFreq1 - $MinorFreq2) ** 2;
  $FracNum = (($SampleSize1 + $SampleSize2)/2) * $SharedNum;
  $FracDen = $SampleSize1 * $SampleSize2 * ((($SampleSize1 + $SampleSize2)/2) - 1);
  $frac = $FracNum / $FracDen;
  $WholeNum = $NumA - $frac;
  $DenFracNum = ($SampleSize1 * $SampleSize2 - ($SampleSize1 + $SampleSize2)/2) * $SharedNum;
  $DenFrac = $DenFracNum / $FracDen;
  $WholeDen = $NumA + $DenFrac;
  $Num12Sum += $WholeNum;
  $Den12Sum += $WholeDen;
  if ($WholeDen != 0){
    $FST = $WholeNum / $WholeDen;
  }
  else{
    $FST = 0;
  }
  if ($FST == 1){
    $FST = $MaxFST;
  }
  push @FSTs, $FST;

#Reynolds FST calculations (pops 1 & 3)
  $SharedNum = $SampleSize1 * ($MinorFreq1 - ($MinorFreq1 ** 2)) + $SampleSize3 * ($MinorFreq3 - ($MinorFreq3 ** 2));
  $NumA = ($MinorFreq1 - $MinorFreq3) ** 2;
  $FracNum = (($SampleSize1 + $SampleSize3)/2) * $SharedNum;
  $FracDen = $SampleSize1 * $SampleSize3 * ((($SampleSize1 + $SampleSize3)/2) - 1);
  $frac = $FracNum / $FracDen;
  $WholeNum = $NumA - $frac;
  $DenFracNum = ($SampleSize1 * $SampleSize3 - ($SampleSize1 + $SampleSize3)/2) * $SharedNum;
  $DenFrac = $DenFracNum / $FracDen;
  $WholeDen = $NumA + $DenFrac;
  $Num13Sum += $WholeNum;
  $Den13Sum += $WholeDen;
  if ($WholeDen != 0){
    $FST = $WholeNum / $WholeDen;
  }
  else{
    $FST = 0;
  }
  if ($FST == 1){
    $FST = $MaxFST;
  }
  if ($WholeDen != 0){
    $FST = $WholeNum / $WholeDen;
  }
  else{
    $FST = 0;
  }
  if ($FST == 1){
    $FST = $MaxFST;
  }
  push @FSTs, $FST;
  
#Reynolds FST calculations (pops 2 & 3)
  $SharedNum = $SampleSize3 * ($MinorFreq3 - ($MinorFreq3 ** 2)) + $SampleSize2 * ($MinorFreq2 - ($MinorFreq2 ** 2));
  $NumA = ($MinorFreq3 - $MinorFreq2) ** 2;
  $FracNum = (($SampleSize3 + $SampleSize2)/2) * $SharedNum;
  $FracDen = $SampleSize3 * $SampleSize2 * ((($SampleSize3 + $SampleSize2)/2) - 1);
  $frac = $FracNum / $FracDen;
  $WholeNum = $NumA - $frac;
  $DenFracNum = ($SampleSize3 * $SampleSize2 - ($SampleSize3 + $SampleSize2)/2) * $SharedNum;
  $DenFrac = $DenFracNum / $FracDen;
  $WholeDen = $NumA + $DenFrac;
  $Num23Sum += $WholeNum;
  $Den23Sum += $WholeDen;
  if ($WholeDen != 0){
    $FST = $WholeNum / $WholeDen;
  }
  else{
    $FST = 0;
  }
  if ($FST == 1){
    $FST = $MaxFST;
  }
  push @FSTs, $FST;
  
#if site passes frequency threshold, record info for site output
  next if ((($MinorFreq1 + $MinorFreq2 + $MinorFreq3) / 3) < $FreqThresh);
  $pos = $site + 1;
  push @positions, $pos;
  push @SiteWins, $window;
  push @SiteFSTAoA, [ @FSTs ];
  push @NucCountAoA, [ @AllNucCounts ];
  $SitePBS = ((-1 * log (1 - $FSTs[0])) + (-1 * log (1 - $FSTs[1])) - (-1 * log (1 - $FSTs[2]))) / 2;
  if ($SitePBS > $MaxPBS){
    $MaxPBS = $SitePBS;
  }
  push @SitePBSs, $SitePBS;
}

#add window data from the last window
push @WinSiteCounts, $SitesAnalyzed;
if ($SitesAnalyzed < $MinCov){
  @pis = (-999, -999, -999);
  push @PiAoA, [ @pis ];
  push @DxyAoA, [ @pis ];
}
else{
  $Pi1Sum = $Pi1Sum / $SitesAnalyzed;
  $Pi2Sum = $Pi2Sum / $SitesAnalyzed;
  $Pi3Sum = $Pi3Sum / $SitesAnalyzed;
  @pis = ($Pi1Sum, $Pi2Sum, $Pi3Sum);
  push @PiAoA, [ @pis ];
  $Dxy12Sum = $Dxy12Sum / $SitesAnalyzed;
  $Dxy13Sum = $Dxy13Sum / $SitesAnalyzed;
  $Dxy23Sum = $Dxy23Sum / $SitesAnalyzed;
  @Dxys = ($Dxy12Sum, $Dxy13Sum, $Dxy23Sum);
  push @DxyAoA, [ @Dxys ];
}
$SitesAnalyzed = 0;
push @MaxPBSs, $MaxPBS;
$MaxPBS = 0;
if (($Den12Sum == 0) || ($Den13Sum == 0) || ($Den23Sum == 0)){
  push @WinPBSs, -999;
  push @WinZs, -999;
  @FSTs = (-999, -999, -999);
  push @FSTAoA, [ @FSTs ];
}
else{
  @FSTs = ();
  $FST = $Num12Sum / $Den12Sum;
  if ($FST >= 1){
    $FST = $MaxFST;
  }
  push @FSTs, $FST;
  $FST = $Num13Sum / $Den13Sum;
  if ($FST >= 1){
    $FST = $MaxFST;
  }
  push @FSTs, $FST;
  $FST = $Num23Sum / $Den23Sum;
  if ($FST >= 1){
    $FST = $MaxFST;
  }
  push @FSTs, $FST;
  push @FSTAoA, [ @FSTs ];
  $WinPBS = ((-1 * log (1 - $FSTs[0])) + (-1 * log (1 - $FSTs[1])) - (-1 * log (1 - $FSTs[2]))) / 2;
  push @WinPBSs, $WinPBS;
  $Z = -1 * log (1 - $FSTs[2]);
  push @WinZs, $Z;
}
print "got data for $chr, window $window\n";

close A;
close B;
close C;


#Quantiles for WinPBS and MaxPBS, plus window medians for PBS and Z
@SortedValues = @WinZs;
@SortedValues = sort { $b <=> $a } @SortedValues;
for ($i = (@SortedValues - 1); $i >= 0; $i--){
  if ($SortedValues[$i] < -0.9){
    splice @SortedValues, $i, 1;
  }
  else{
    last;
  }
}
if ((@SortedValues % 2) == 0){
  $i = @SortedValues / 2;
  $MedianZ = ($SortedValues[$i] + $SortedValues[$i-1]) / 2;
}
else{
  $i = (@SortedValues / 2) - 1;
  $MedianZ = $SortedValues[$i];
}

@SortedValues = @WinPBSs;
@SortedValues = sort { $b <=> $a } @SortedValues;
for ($i = (@SortedValues - 1); $i >= 0; $i--){
  if ($SortedValues[$i] == -999){
    splice @SortedValues, $i, 1;
  }
  else{
    last;
  }
}
for ($i = 0; $i < @WinPBSs; $i++){
  if ($WinPBSs[$i] == -999){
    push @WinPBSQuantiles, -999;
    next;
  }
  for ($j = 0; $j < @SortedValues; $j++){
    if ($WinPBSs[$i] >= $SortedValues[$j]){
      $quantile = ($j / @SortedValues);
      push @WinPBSQuantiles, $quantile;
      last;
    }
  }
}

if ((@SortedValues % 2) == 0){
  $i = @SortedValues / 2;
  $MedianPBS = ($SortedValues[$i] + $SortedValues[$i-1]) / 2;
}
else{
  $i = (@SortedValues / 2) - 1;
  $MedianPBS = $SortedValues[$i];
}

@SortedValues = @MaxPBSs;
@SortedValues = sort { $b <=> $a } @SortedValues;
for ($i = (@SortedValues - 1); $i >= 0; $i--){
  if ($SortedValues[$i] == -999){
    splice @SortedValues, $i, 1;
  }
  else{
    last;
  }
}
for ($i = 0; $i < @MaxPBSs; $i++){
  if ($MaxPBSs[$i] == -999){
    push @MaxPBSQuantiles, -999;
    next;
  }
  for ($j = 0; $j < @SortedValues; $j++){
    if ($MaxPBSs[$i] >= $SortedValues[$j]){
      $quantile = ($j / @SortedValues);
      push @MaxPBSQuantiles, $quantile;
      last;
    }
  }
}

#calculate window PBEs
for ($i = 0; $i < @WinPBSs; $i++){
  if ($WinPBSs[$i] == -999){
    push @WinPBEs, -999;
  }
  else{
    $WinPBE = $WinPBSs[$i] - ($WinZs[$i] * ($MedianPBS / $MedianZ));
    push @WinPBEs, $WinPBE;
  }
}

#calculate site PBEs
for ($i = 0; $i < @SitePBSs; $i++){
  $Z = -1 * log (1 - $SiteFSTAoA[$i][2]);
  $SitePBE = $SitePBSs[$i] - ($Z * ($MedianPBS / $MedianZ));
  push @SitePBEs, $SitePBE;
}

#get max PBE for each window
$j = 0;
for ($i = 0; $i < @WinPBEs; $i++){
  $MaxPBE = -999;
  while (($j < @SiteWins) && ($SiteWins[$j] == $i)){
    if ($SitePBEs[$j] > $MaxPBE){
      $MaxPBE = $SitePBEs[$j];
    }
    $j++;
  }
  push @MaxPBEs, $MaxPBE;
}

#add quantiles for WinPBE and MaxPBE
@SortedValues = @WinPBEs;
@SortedValues = sort { $b <=> $a } @SortedValues;
for ($i = (@SortedValues - 1); $i >= 0; $i--){
  if ($SortedValues[$i] == -999){
    splice @SortedValues, $i, 1;
  }
  else{
    last;
  }
}
for ($i = 0; $i < @WinPBEs; $i++){
  if ($WinPBEs[$i] == -999){
    push @WinPBEQuantiles, -999;
    next;
  }
  for ($j = 0; $j < @SortedValues; $j++){
    if ($WinPBEs[$i] >= $SortedValues[$j]){
      $quantile = ($j / @SortedValues);
      push @WinPBEQuantiles, $quantile;
      last;
    }
  }
}
@SortedValues = @MaxPBEs;
@SortedValues = sort { $b <=> $a } @SortedValues;
for ($i = (@SortedValues - 1); $i >= 0; $i--){
  if ($SortedValues[$i] == -999){
    splice @SortedValues, $i, 1;
  }
  else{
    last;
  }
}
for ($i = 0; $i < @MaxPBEs; $i++){
  if ($MaxPBEs[$i] == -999){
    push @MaxPBEQuantiles, -999;
    next;
  }
  for ($j = 0; $j < @SortedValues; $j++){
    if ($MaxPBEs[$i] >= $SortedValues[$j]){
      $quantile = ($j / @SortedValues);
      push @MaxPBEQuantiles, $quantile;
      last;
    }
  }
}

#Output
open O, ">$OutputFile" or die;
print O "Chr\tWinStart\tWinStop\tSitesAnalyzed\tWinPBE\tWinPBEQuantile\tMaxPBEs\tMaxPBEQuantiles\tWinPBSs\tWinPBSQuantiles\tMaxPBSs\tMaxPBSQuantiles\tpi_$pops[0]\tpi_$pops[1]\tpi_$pops[2]\tDxy_$pops[0]_$pops[1]\tDxy_$pops[0]_$pops[2]\tDxy_$pops[1]_$pops[2]\tFST_$pops[0]_$pops[1]\tFST_$pops[0]_$pops[2]\tFST_$pops[1]_$pops[2]\tMedianPBS:\t$MedianPBS\tMedianZ:\t$MedianZ\n";
 for ($i = 0; $i < @WinPBEs; $i++){
  print O "$chr\t$WinStarts[$i]\t$WinStops[$i]\t$WinSiteCounts[$i]\t$WinPBEs[$i]\t$WinPBEQuantiles[$i]\t$MaxPBEs[$i]\t$MaxPBEQuantiles[$i]\t$WinPBSs[$i]\t$WinPBSQuantiles[$i]\t$MaxPBSs[$i]\t$MaxPBSQuantiles[$i]\t$PiAoA[$i][0]\t$PiAoA[$i][1]\t$PiAoA[$i][2]\t$DxyAoA[$i][0]\t$DxyAoA[$i][1]\t$DxyAoA[$i][2]\t$FSTAoA[$i][0]\t$FSTAoA[$i][1]\t$FSTAoA[$i][2]\n";
}
close O;

open S, ">$SNPFile" or die;
print S "SNP_pos\tPBE\tPBS\tFST_$pops[0]_$pops[1]\tFST_$pops[0]_$pops[2]\tFST_$pops[1]_$pops[2]\t$pops[0]_As\t$pops[0]_Cs\t$pops[0]_Gs\t$pops[0]_Ts\t$pops[1]_As\t$pops[1]_Cs\t$pops[1]_Gs\t$pops[1]_Ts\t$pops[2]_As\t$pops[2]_Cs\t$pops[2]_Gs\t$pops[2]_Ts\tMedianPBS:\t$MedianPBS\tMedianZ:\t$MedianZ\n";
for ($i = 0; $i < @SitePBSs; $i++){
  print S "$positions[$i]\t$SitePBEs[$i]\t$SitePBSs[$i]\t$SiteFSTAoA[$i][0]\t$SiteFSTAoA[$i][1]\t$SiteFSTAoA[$i][2]\t$NucCountAoA[$i][0]\t$NucCountAoA[$i][1]\t$NucCountAoA[$i][2]\t$NucCountAoA[$i][3]\t$NucCountAoA[$i][4]\t$NucCountAoA[$i][5]\t$NucCountAoA[$i][6]\t$NucCountAoA[$i][7]\t$NucCountAoA[$i][8]\t$NucCountAoA[$i][9]\t$NucCountAoA[$i][10]\t$NucCountAoA[$i][11]\n";
}
close S;
