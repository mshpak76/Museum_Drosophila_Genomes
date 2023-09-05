#!/usr/bin/perl -w
use strict;

#Identifies outlier windows for a previously calculated window statistic, groups nearby outlier windows into outlier regions, and identifies nearest gene, etc.
#This version excludes windows with statistic value "NA", and does not permit outlier regions to span across missing windows

my $InputHeaderRows = 1;  #how many header rows are present in the input file before the data
my $InputColumn = 4;  #which column from input file is the statistic in, starting from column 0
#assumes that column 0 is chromosome arm, columns 1 and 2 are window start and stop positions. Column 8 is for PBS

my $OutlierQuantile = 0.01; #the proportion of windows on each chromosome defined as outliers, set to 0.01, try 0.05 for comparison
my $OutliersHigh = 1;  #set to 1 if looking for high outliers, -1 if looking for low

my $MaxBtwnOutliers = 5;  #maximum number of non-outlier windows between two outlier windows that will be grouped into one outlier region  (1 means only use neighboring windows)



#my @InputFiles = ('Nineteenth_ModernLund_TwentiethZI_WinPBE_ChrX.txt','Nineteenth_ModernLund_TwentiethZI_WinPBE_Chr2L.txt','Nineteenth_ModernLund_TwentiethZI_WinPBE_Chr2R.txt','Nineteenth_ModernLund_TwentiethZI_WinPBE_Chr3L.txt', 'Nineteenth_ModernLund_TwentiethZI_WinPBE_Chr3R.txt');

#my @InputFiles = ('Nineteenth_ModernLund_FRZI_WinPBE_ChrX.txt','Nineteenth_ModernLund_FRZI_WinPBE_Chr2L.txt','Nineteenth_ModernLund_FRZI_WinPBE_Chr2R.txt','Nineteenth_ModernLund_FRZI_WinPBE_Chr3L.txt', 'Nineteenth_ModernLund_FRZI_WinPBE_Chr3R.txt');

#my @InputFiles =('Nineteenth_ModernLund_StockholmZI_WinPBE_ChrX.txt','Nineteenth_ModernLund_StockholmZI_WinPBE_Chr2L.txt','Nineteenth_ModernLund_StockholmZI_WinPBE_Chr2R.txt','Nineteenth_ModernLund_StockholmZI_WinPBE_Chr3L.txt','Nineteenth_ModernLund_StockholmZI_WinPBE_Chr3R.txt');

#my @InputFiles = ('Nineteenth_FR_StockholmZI_WinPBE_ChrX.txt','Nineteenth_FR_StockholmZI_WinPBE_Chr2L.txt','Nineteenth_FR_StockholmZI_WinPBE_Chr2R.txt','Nineteenth_FR_StockholmZI_WinPBE_Chr3L.txt','Nineteenth_FR_StockholmZI_WinPBE_Chr3R.txt');


my @InputFiles = ('Nineteenth_ModernLund_LyonZI_WinPBE_Chr2L.txt','Nineteenth_ModernLund_LyonZI_WinPBE_Chr2R.txt','Nineteenth_ModernLund_LyonZI_WinPBE_Chr3L.txt','Nineteenth_ModernLund_LyonZI_WinPBE_Chr3R.txt', 'Nineteenth_ModernLund_LyonZI_WinPBE_ChrX.txt');

#my @InputFiles = ('Twentieth_ModernLund_LyonZI_WinPBE_ChrX.txt','Twentieth_ModernLund_LyonZI_WinPBE_Chr2L.txt','Twentieth_ModernLund_LyonZI_WinPBE_Chr2R.txt','Twentieth_ModernLund_LyonZI_WinPBE_Chr3L.txt', 'Twentieth_ModernLund_LyonZI_WinPBE_Chr3R.txt');

#my @InputFiles = ('Museum_ModernLund_FRZI_WinPBE_ChrX.txt','Museum_ModernLund_FRZI_WinPBE_Chr2L.txt','Museum_ModernLund_FRZI_WinPBE_Chr2R.txt','Museum_ModernLund_FRZI_WinPBE_Chr3L.txt', 'Museum_ModernLund_FRZI_WinPBE_Chr3R.txt');

#my @InputFiles = ('Twentieth_Stockholm_FRZI_WinPBE_ChrX.txt','Twentieth_Stockholm_FRZI_WinPBE_Chr2L.txt','Twentieth_Stockholm_FRZI_WinPBE_Chr2R.txt','Twentieth_Stockholm_FRZI_WinPBE_Chr3L.txt','Twentieth_Stockholm_FRZI_WinPBE_Chr3R.txt');

#my @InputFiles =('Twentieth_Stockholm_ModernLundZI_WinPBE_Chr2L.txt','Twentieth_Stockholm_ModernLundZI_WinPBE_Chr2R.txt','Twentieth_Stockholm_ModernLundZI_WinPBE_Chr3L.txt','Twentieth_Stockholm_ModernLundZI_WinPBE_Chr3R.txt','Twentieth_Stockholm_ModernLundZI_WinPBE_ChrX.txt');

#my @InputFiles =('Museum_Stockholm_ModernLundZI_WinPBE_Chr2L.txt','Museum_Stockholm_ModernLundZI_WinPBE_Chr2R.txt','Museum_Stockholm_ModernLundZI_WinPBE_Chr3L.txt','Museum_Stockholm_ModernLundZI_WinPBE_Chr3R.txt','Museum_Stockholm_ModernLundZI_WinPBE_ChrX.txt');

#my @InputFiles =('Museum_Stockholm_FRZI_WinPBE_Chr2L.txt','Museum_Stockholm_FRZI_WinPBE_Chr2R.txt','Museum_Stockholm_FRZI_WinPBE_Chr3L.txt','Museum_Stockholm_FRZI_WinPBE_Chr3R.txt','Museum_Stockholm_FRZI_WinPBE_ChrX.txt');

#my @InputFiles = ('Nineteenth_LundOnly_ModernLund_TwentiethZI_WinPBE_Chr2L.txt','Nineteenth_LundOnly_ModernLund_TwentiethZI_WinPBE_Chr2R.txt','Nineteenth_LundOnly_ModernLund_TwentiethZI_WinPBE_Chr3R.txt','Nineteenth_LundOnly_ModernLund_TwentiethZI_WinPBE_Chr3L.txt','Nineteenth_LundOnly_ModernLund_TwentiethZI_WinPBE_ChrX.txt');

#my @InputFiles = ('Nineteenth_LundOnly_ModernLund_FRZI_WinPBE_Chr2L.txt','Nineteenth_LundOnly_ModernLund_FRZI_WinPBE_Chr2R.txt','Nineteenth_LundOnly_ModernLund_FRZI_WinPBE_Chr3L.txt','Nineteenth_LundOnly_ModernLund_FRZI_WinPBE_Chr3R.txt','Nineteenth_LundOnly_ModernLund_FRZI_WinPBE_ChrX.txt');

#my @InputFiles = ('Nineteenth_LundOnly_ModernLund_StockholmZI_WinPBE_Chr2L.txt','Nineteenth_LundOnly_ModernLund_StockholmZI_WinPBE_Chr2R.txt', 'Nineteenth_LundOnly_ModernLund_StockholmZI_WinPBE_Chr3L.txt','Nineteenth_LundOnly_ModernLund_StockholmZI_WinPBE_Chr3R.txt','Nineteenth_LundOnly_ModernLund_StockholmZI_WinPBE_Chr2R.txt', 'Nineteenth_LundOnly_ModernLund_StockholmZI_WinPBE_ChrX.txt');

#my $OutputFile = 'Nineteenth_ModernLund_Twentieth_PBE_regions_01.txt';
#my $OutputFile = 'NineteenthCent_LundFR_PBE_regions_01.txt';

#my $OutputFile = 'NineteenthCent_LundStockholm_PBE_regions_01.txt';
#my $OutputFile = 'NineteenthCent_FRStockholm_PBE_regions_01.txt';
my $OutputFile = 'NineteentCent_Lund_Lyon_PBE_regions_01.txt';

#my $OutputFile = 'TwentiethCent_Lund_Lyon_PBE_regions_01.txt';
#my $OutputFile = 'Museum_LundFR_PBE_regions_01.txt';
#my $OutputFile = 'TwentiethCent_StockholmFR_PBE_regions_01.txt';
#my $OutputFile = 'TwentiethCent_StockholmModernLund_PBE_regions_01.txt';
#my $OutputFile = 'Museum_StockholmLund_PBE_regions_01.txt';
#my $OutputFile = 'Museum_StockholmFR_PBE_regions_01.txt';

#my $OutputFile = 'NineteenthLundOnly_LundTwentiethPBE_regions_01.txt';
#my $OutputFile = 'NineteenthLundOnly_FRPBE_regions_01.txt';
#my $OutputFile = 'NineteenthLundOnly_StockholmPBE_regions_01.txt';

my $ExonFilePrefix = 'annot550_'; #assumes full name of annotation file is something like annot_chr2L.txt

my $MissingSymbol = '-999';

#declaring variables (no need to modify these)
my $f = 0;
my $file = '';
my $i = 0;
my $j = 0;
my $g = 0;
my $ThreshValue = 0;
my $chr = '';
my $LastSig = 0;
my $RegionStart = 0;
my $RegionStop = 0;
my $RegionSigs = 0;
my $RegionMax = 0;
my $RegionMaxCtr = 0;
my $FirstSig = 0;
my $RegionWins = 0;

my @line = ();
my @chr = ();
my @WinStart = ();
my @WinStop = ();
my @statistic = ();
my @SortedValues = ();
my @RegionChr = ();
my @RegionStart = ();
my @RegionStop = ();
my @RegionMax = ();
my @RegionMaxCtr = ();
my @RegionSigs = ();
my @RegionWins = ();

my @OutputAoA = ();
my @CDSAoA = ();
my @ClosestGeneList = ();
my @ClosestFbgnList = ();
my @GeneOrientationList = ();
my @GenePositionList = ();
my $k = 0;
my $h = 0;
#my $NearestExon = -1;
my $NearestDistance = 1000000000;
my $distance = 0;
#my @GenesWithin = ();
my @UpstreamGenes = ();
my @DownstreamGenes = ();
my @ClosestGenes = ();
my $AlreadyListed = 0;
my $orientation = "";
my $GenePosition = 0;
my $UpstreamLimit = 0;
my $DownstreamLimit = 0;
my $direction = 0;
my $GeneList = "";
my @ClosestFbgns = ();
my $FbgnList = "";
my @OverlapGenes = ();
my @OverlapFbgns = ();
my @OverlapGeneList = ();
my @OverlapFbgnList = ();
my $ClosestGene = "";
my $ClosestDist = 1000000000;

#Get input for statistic in windows
for ($f = 0; $f < @InputFiles; $f++){
  $file = $InputFiles[$f];
  open I, "<$file" or die "can not open $file\n";
  for ($i = 0; $i < $InputHeaderRows; $i++){
    scalar (<I>);
  }
  @chr = ();
  @WinStart = ();
  @WinStop = ();
  @statistic = ();
  while (<I>){
    chomp;
    last if m/^$/;
    @line = split;
    next if ($line[$InputColumn] =~ $MissingSymbol);
    push @chr, $line[0];
    $chr = $line[0];
    push @WinStart, $line[1];
    push @WinStop, $line[2];
    push @statistic, $line[$InputColumn];
  }
  close I;

#Calculate threshold values
  if ($OutliersHigh == 1){
    @SortedValues = sort{$b<=>$a} @statistic;
  }
  else{
    @SortedValues = sort{$a<=>$b} @statistic;
  }
  $ThreshValue = int(@SortedValues * $OutlierQuantile - 1);
  $ThreshValue = $SortedValues[$ThreshValue];
  print "For $chr[0], calculated threshold value of $ThreshValue\n";

#Define "outlier regions" 

  my $thresh = 0;
  
  for ($i = 0; $i < @statistic; $i++){


    if ( ($statistic[$i] * $OutliersHigh) > ($ThreshValue * $OutliersHigh) ){
      $FirstSig = $i;
      $RegionSigs = 1;
      $RegionMax = $statistic[$i];
      $RegionMaxCtr = ($WinStart[$i] + $WinStop[$i]) / 2;
      if ($chr ne $chr[$i] ){
	die "Found different chromosome arm names in same input file ($chr vs $chr[$i]).  Use only one input file per chromosome arm.\n";
      }
      $chr = $chr[$i];
      $RegionStart = $WinStart[$i];
      $RegionStop = $WinStop[$i];
      $LastSig = $i;
      $i++;
      while ( $i <= ($LastSig + $MaxBtwnOutliers) ){  
	if ( ($statistic[$i] * $OutliersHigh) > ($ThreshValue * $OutliersHigh) ){
	  $RegionSigs++;
	  $RegionStop = $WinStop[$i];
	  $LastSig = $i;
	  if ( ($statistic[$i] * $OutliersHigh) > ($RegionMax * $OutliersHigh) ){
	    $RegionMax = $statistic[$i];
	    $RegionMaxCtr = ($WinStart[$i] + $WinStop[$i]) / 2;
	  }
	}
	last if (($i == (@WinStop - 1)) || (($WinStop[$i] + 1) != $WinStart[$i+1]));
	$i++;
	last if ($i >= @statistic);
	#	last if ($chr[$i] ne $chr);
      }
      $i = $LastSig;
      $RegionWins = ($LastSig - $FirstSig) + 1;
      push @RegionChr, $chr;
      push @RegionStart, $RegionStart;
      push @RegionStop, $RegionStop;
      push @RegionMax, $RegionMax;
      push @RegionMaxCtr, $RegionMaxCtr;
      push @RegionSigs, $RegionSigs;
      push @RegionWins, $RegionWins;
    }
  }
  $i = @RegionStart;
  print "Total outlier regions on $chr: $i\n";

#get exon data for this arm
  @CDSAoA = ();
  $file = $ExonFilePrefix . $chr . '.txt';
  open E, "<$file" or die "can not open exon file $file\n";
  while (<E>){
    chomp;
    last if m/^$/;
    @line = split;
    push @CDSAoA, [ @line ];
  }
  close E;

#Begin the process of identifying closest gene to outlier region center
  $j = 0;
  for ($i = 0; $i < @RegionMaxCtr; $i++){
#  $NearestExon = -1;
    $NearestDistance = 1000000000;
    $orientation = 'n';
    $GenePosition = 9000000000;
    #  @GenesWithin = ();
    @UpstreamGenes = ();
    @DownstreamGenes = ();
    @ClosestGenes = ();
    @ClosestFbgns = ();
    $UpstreamLimit = 1000000000;
    $DownstreamLimit = 0;
    $GeneList = "";
    @OverlapFbgns = ();
    @OverlapGenes = ();
    
    print "Evaluating outlier region $i.  $RegionChr[$i]:$RegionMaxCtr[$i]\n";

#Move $j to the first exon that starts after the sweep center position
    while ( $RegionMaxCtr[$i] > $CDSAoA[$j][4] ){
      $j++;
      last if ($j == (@CDSAoA - 1));
    }   
    while ( $RegionMaxCtr[$i] < $CDSAoA[$j][4] ){
      last if ($j == 0);
      $j--;
    }
    $k = $j;
    
#Scan upstream 500 exons to check if sweep center is within any exons, or which exon(s) are closest
    while ($j > ($k - 500) ){
      
#Check if we're inside an exon
      if ($RegionMaxCtr[$i] < $CDSAoA[$j][5]){
	if ($CDSAoA[$j][6] eq 'e'){
	  $orientation = 'e';
	}
	elsif ($orientation eq 'n'){
	  $orientation = $CDSAoA[$j][6];
	}
	if ($NearestDistance == 0){
	  $AlreadyListed = 0;
	  for ($g = 0; $g < @ClosestGenes; $g++){
	    if ($ClosestGenes[$g] eq $CDSAoA[$j][1]){
	      $AlreadyListed = 1;
	      last;
	    }
	  }
	  if ($AlreadyListed == 0){
	    if ($CDSAoA[$j][6] eq 'e'){
	      unshift @ClosestGenes, $CDSAoA[$j][1];
	      unshift @ClosestFbgns, $CDSAoA[$j][0];
	    }
	    else{
	      push @ClosestGenes, $CDSAoA[$j][1];
	      push @ClosestFbgns, $CDSAoA[$j][0];
	    }
	  }
	}
	else{
	  @ClosestGenes = ();
	  @ClosestFbgns = ();
	  push @ClosestGenes, $CDSAoA[$j][1];
	  push @ClosestFbgns, $CDSAoA[$j][0];
	  $NearestDistance = 0;
	}
      }

#Otherwise, check if this is closest exon so far, and add to list of upstream genes
      else{
	$distance = $RegionMaxCtr[$i] - $CDSAoA[$j][5];
	if ($distance < $NearestDistance){
	  @ClosestGenes = ();
	  @ClosestFbgns = ();
	  push @ClosestGenes, $CDSAoA[$j][1];
	  push @ClosestFbgns, $CDSAoA[$j][0];
	  $NearestDistance = $distance;
	}
	elsif ($distance == $NearestDistance){
	  $AlreadyListed = 0;
	  for ($g = 0; $g < @ClosestGenes; $g++){
	    if ($ClosestGenes[$g] eq $CDSAoA[$j][1]){
	      $AlreadyListed = 1;
	      last;
	    }
	  }
	  if ($AlreadyListed == 0){
	    push @ClosestGenes, $CDSAoA[$j][1];
	    push @ClosestFbgns, $CDSAoA[$j][0];
	  }
	}
	
	$AlreadyListed = 0;
	for ($g = 0; $g < @UpstreamGenes; $g++){
	  if ($UpstreamGenes[$g] eq $CDSAoA[$j][1]){
	   $AlreadyListed = 1;
	   last;
	 }
	}
	if ($AlreadyListed == 0){
	  push @UpstreamGenes, $CDSAoA[$j][1];
	}
      }
      
      $j--;
      last if ($j < 0);
    }

#Scan downstream 500 exons and check which exon(s) are closest
    $j = $k + 1;
    while ($j < ($k + 500) ){
      $distance = $CDSAoA[$j][4] - $RegionMaxCtr[$i];
      if ($distance < $NearestDistance){
	@ClosestGenes = ();
	@ClosestFbgns = ();
	push @ClosestGenes, $CDSAoA[$j][1];
	push @ClosestFbgns, $CDSAoA[$j][0];
	$NearestDistance = $distance;
      }
      elsif ($distance == $NearestDistance){
	$AlreadyListed = 0;
	for ($g = 0; $g < @ClosestGenes; $g++){
	  if ($ClosestGenes[$g] eq $CDSAoA[$j][1]){
	    $AlreadyListed = 1;
	    last;
	  }
	}
	if ($AlreadyListed == 0){
	  push @ClosestGenes, $CDSAoA[$j][1];
	  push @ClosestFbgns, $CDSAoA[$j][0];
	}
      }
      
      $AlreadyListed = 0;
      for ($g = 0; $g < @DownstreamGenes; $g++){
	if ($DownstreamGenes[$g] eq $CDSAoA[$j][1]){
	  $AlreadyListed = 1;
	  last;
	}
      }
      if ($AlreadyListed == 0){
	push @DownstreamGenes, $CDSAoA[$j][1];
      }
      $j++;
      last if ($j == @CDSAoA);
    }

#Determine the boundaries of the closest gene (scan 500 exons up and down from start position)
    $UpstreamLimit = 1000000000;
    $DownstreamLimit = 0;
    for ($j = $k - 500; $j < $k + 500; $j++){
      next if ($j < 0);
      if ($CDSAoA[$j][1] eq $ClosestGenes[0]){
	if ($CDSAoA[$j][4] < $UpstreamLimit){
	  $direction = $CDSAoA[$j][3];
	  $UpstreamLimit = $CDSAoA[$j][4];
	}
	if ($CDSAoA[$j][5] > $DownstreamLimit){
	  $DownstreamLimit = $CDSAoA[$j][5];
	}
      }
      last if ( ($j + 1) == @CDSAoA );
    }

#If within gene, resolve exon vs. intron, and record gene position on interval 0 to 1 (5' to 3')
    if ( ($RegionMaxCtr[$i] > $UpstreamLimit) && ($RegionMaxCtr[$i] < $DownstreamLimit) ){
      if ($orientation eq 'n'){
	$orientation = 'i';
      }
      $GenePosition = ($RegionMaxCtr[$i] - $UpstreamLimit) / ($DownstreamLimit - $UpstreamLimit);
      if ($direction ne '+'){
	$GenePosition = 1 - $GenePosition;
      }
    }

#If outside gene, note up or downstream (accounting for gene's strand), assign - or + GenePosition (bp away) accordingly
    elsif ($RegionMaxCtr[$i] < $UpstreamLimit){
      if ($direction eq '+'){
	$orientation = 'u';
	$GenePosition = $RegionMaxCtr[$i] - $UpstreamLimit;
      }
      else{
	$orientation = 'd';
	$GenePosition = $UpstreamLimit - $RegionMaxCtr[$i];
      }
    }

    else{
      if ($direction eq '+'){
	$orientation = 'd';
	$GenePosition = $RegionMaxCtr[$i] - $DownstreamLimit;
      }
      else{
	$orientation = 'u';
	$GenePosition = $DownstreamLimit - $RegionMaxCtr[$i];
      }
    }

#Obtain all genes that overlap the full outlier region, plus the next exon to each side
    for ($j = 0; $j < @CDSAoA; $j++){
      next if ($RegionStart[$i] > $CDSAoA[$j][4]);
#if we hit our first exon that starts after the region start, go back and get genes with exons that overlap this boundary, or else the closest exon before it
      if (@OverlapFbgns == 0){
	$ClosestGene = '';
	$ClosestDist = 1000000000;
	for ($k = $j - 1; $k >= 0; $k--){
	  if ($RegionStart[$i] < $CDSAoA[$k][5]){
	    for ($g = 0; $g < @OverlapFbgns; $g++){
	      if ($CDSAoA[$k][0] eq $OverlapFbgns[$g]){
		splice @OverlapFbgns, $g, 1;
		splice @OverlapGenes, $g, 1;
		last;
	      }
	    }
	    push @OverlapFbgns, $CDSAoA[$k][0];
	    push @OverlapGenes, $CDSAoA[$k][1];
	    next;
	  }
	  if ( ($RegionStart[$i] - $CDSAoA[$k][5]) < $ClosestDist){
	    $ClosestGene = $k;
	    $ClosestDist = $RegionStart[$i] - $CDSAoA[$k][5];
	  }
	  if ( ($k <= ($j-500)) || ($k == 0) ){
	    if (@OverlapFbgns == 0){
	       push @OverlapFbgns, $CDSAoA[$ClosestGene][0];
	       push @OverlapGenes, $CDSAoA[$ClosestGene][1];
	     }
	    last;
	  }
	}
	$ClosestGene = '';
	$ClosestDist = 1000000000;
      }
#Gather all genes with exons that fall completely within the outlier region
      if ($CDSAoA[$j][5] < $RegionStop[$i]){
	for ($g = 0; $g < @OverlapFbgns; $g++){
	  if ($CDSAoA[$j][0] eq $OverlapFbgns[$g]){
	    splice @OverlapFbgns, $g, 1;
	    splice @OverlapGenes, $g, 1;
	    last;
	  }
	}
	push @OverlapFbgns, $CDSAoA[$j][0];
	push @OverlapGenes, $CDSAoA[$j][1];
	next;
      }
#when we reach the first exon that starts after the region stop, include this gene if there were none overlapping the outlier region stop, then terminate the loop
      elsif ($CDSAoA[$j][4] > $RegionStop[$i]){
	if ($ClosestDist > 0){
	  for ($g = 0; $g < @OverlapFbgns; $g++){
	    if ($CDSAoA[$j][0] eq $OverlapFbgns[$g]){
	      splice @OverlapFbgns, $g, 1;
	      splice @OverlapGenes, $g, 1;
	      last;
	    }
	  }
	  push @OverlapFbgns, $CDSAoA[$j][0];
	  push @OverlapGenes, $CDSAoA[$j][1];
	}
	last;
      }
#include genes with exons that overlap the region stop
      else{
	for ($g = 0; $g < @OverlapFbgns; $g++){
	  if ($CDSAoA[$j][0] eq $OverlapFbgns[$g]){
	    splice @OverlapFbgns, $g, 1;
	    splice @OverlapGenes, $g, 1;
	    last;
	  }
	}
	push @OverlapFbgns, $CDSAoA[$j][0];
	push @OverlapGenes, $CDSAoA[$j][1];
	$ClosestDist = 0;
      }
    }

#Record information
    push @GeneOrientationList, $orientation;
    push @GenePositionList, $GenePosition;
    if (@ClosestGenes > 1){
      $GeneList = $ClosestGenes[0];
      $FbgnList = $ClosestFbgns[0];
      for ($g = 1; $g < @ClosestGenes; $g++){
	$GeneList = $GeneList . '/' . $ClosestGenes[$g];
	$FbgnList = $FbgnList . '/' . $ClosestFbgns[$g];
      }
      push @ClosestGeneList, $GeneList;
      push @ClosestFbgnList, $FbgnList;
    }
    else{
      push @ClosestGeneList, $ClosestGenes[0];
      push @ClosestFbgnList, $ClosestFbgns[0];
    }
    if (@OverlapGenes > 1){
      $GeneList = $OverlapGenes[0];
      $FbgnList = $OverlapFbgns[0];
      for ($g = 1; $g < @OverlapGenes; $g++){
	$GeneList = $GeneList . '/' . $OverlapGenes[$g];
	$FbgnList = $FbgnList . '/' . $OverlapFbgns[$g];
      }
      push @OverlapGeneList, $GeneList;
      push @OverlapFbgnList, $FbgnList;
    }
    elsif (@OverlapGenes == 1){
      push @OverlapGeneList, $OverlapGenes[0];
      push @OverlapFbgnList, $OverlapFbgns[0];
    }
    else{
      push @OverlapGeneList, '';
      push @OverlapFbgnList, '';
    }
  }
  for ($i = 0; $i < @RegionMaxCtr; $i++){
    @line = ();
    push @line, $RegionChr[$i];
    push @line, $RegionStart[$i];
    push @line, $RegionStop[$i];
    push @line, $RegionWins[$i];
    push @line, $RegionSigs[$i];
    push @line, $RegionMax[$i];
    push @line, $RegionMaxCtr[$i];
    push @line, $ClosestFbgnList[$i];
    push @line, $ClosestGeneList[$i];
    push @line, $GeneOrientationList[$i];
    push @line, $GenePositionList[$i];
    push @line, $OverlapFbgnList[$i];
    push @line, $OverlapGeneList[$i];
    push @OutputAoA, [ @line ];
  }
    @RegionChr = ();
    @RegionStart = ();
    @RegionStop = ();
    @RegionMax = ();
    @RegionMaxCtr = ();
    @ClosestGeneList = ();
    @ClosestFbgnList = ();
    @GeneOrientationList = ();
    @GenePositionList = ();
    @OverlapGeneList = ();
    @OverlapFbgnList = ();
  @RegionSigs = ();
  @RegionWins = ();
}

#Output
open O, ">$OutputFile";
print O "chr\tstart_pos\tstop_pos\t#total_windows\t#signif_windows\tMax_Stat\tmax_window_center\tnearest_gene(s)\tnearest_Fbgn(s)\torientation\tposition_vs_gene\tOverlapFbgns\tOverlapGenes\n";
for ($i = 0; $i < @OutputAoA; $i++){
  for ($j = 0; $j < @{$OutputAoA[$i]}; $j++){
    print O $OutputAoA[$i][$j];
    if ($j < (@{$OutputAoA[$i]} - 1) ){
      print O "\t";
    }
    else{
      print O "\n";
    }
  }
}
