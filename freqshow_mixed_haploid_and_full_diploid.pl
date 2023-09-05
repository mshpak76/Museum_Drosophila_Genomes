#!/usr/bin/perl -w
#use strict;

#Create a file that gives allele counts at every site for a single population (A, C, G, T), for one chromosome arm

my $PopLabel = 'H';

my @InLines = ('H4','H5','H6','H12','H10','H11','H13','H9','H25');  #List haploids first
#hap H4,5,6,12
#dip 10,11,13,9,25 for autosome

my $FirstDiploid = 4;  #Haploids start from element 0 above

my @chrs = ('Chr2L','Chr2R','Chr3L','Chr3R');

my @InHandles = @InLines;
for ($i = 0; $i < @InHandles; $i++){
  $InHandles[$i] =~ s/1/A/;
  $InHandles[$i] =~ s/2/B/;
  $InHandles[$i] =~ s/3/C/;
  $InHandles[$i] =~ s/4/D/;
  $InHandles[$i] =~ s/5/E/;
  $InHandles[$i] =~ s/6/F/;
  $InHandles[$i] =~ s/7/G/;
  $InHandles[$i] =~ s/8/H/;
  $InHandles[$i] =~ s/9/I/;
  $InHandles[$i] =~ s/0/J/;
}

#declarations
my $c = 0;
my $f = 0;
my $i = 0;
my $j = 0;
my $s = 0;
my $chr = '';
my $kb = -1;
my $file = '';
my $file1 = '';
my $path = '';
my $CountFile = '';
my @ChrLines = ();
my @ChrHandles = ();
my @line = ();
my @InputAoA = ();
my @SiteCounts = ();

#Chr loop (not indented)
for ($c = 0; $c < @chrs; $c++){
$chr = $chrs[$c];
$kb = -1;
$CountFile = 'AllCounts_' . $PopLabel . '_NoInv_' . $chr . '.txt';
open O, ">$CountFile" or die;

#ind loop - open files
@ChrLines = ();
@ChrHandles = ();
for ($f = 0; $f < @InLines; $f++){
  $file = $InLines[$f] . "_" . $chr . "\.fas1k";
  $path = "$file";
  if (-e $path){
    open $InHandles[$f], "$file";
    push @ChrLines, $InLines[$f];
    push @ChrHandles, $InHandles[$f];
  } 
}
$file1 = $ChrHandles[0];
while (<$file1>){
  @InputAoA = ();
  chomp;
  $kb++;
  last if m/^$/;
  @line = split(//, $_);
  push @InputAoA, [ @line ];
  for ($f = 1; $f < @ChrLines; $f++){
    $file = $ChrHandles[$f];
    $_ = (<$file>);
    chomp;
    @line = split(//, $_);
    push @InputAoA, [ @line ];
  }
  if (($kb % 1000) == 0){
    print "kb $kb\n";
  }

#Obtain and record allele counts for each site  
  for ($s = 0; $s < @{$InputAoA[0]}; $s++){
    @SiteCounts = (0,0,0,0);
    for ($i = 0; $i < $FirstDiploid; $i++){
      if ($InputAoA[$i][$s] eq 'A'){
	$SiteCounts[0]++;
      }
      elsif ($InputAoA[$i][$s] eq 'C'){
	$SiteCounts[1]++;
      }
      elsif ($InputAoA[$i][$s] eq 'G'){
	$SiteCounts[2]++;
      }
      elsif ($InputAoA[$i][$s] eq 'T'){
	$SiteCounts[3]++;
      }
    }
    for ($i = $FirstDiploid; $i < @InputAoA; $i++){
      if ($InputAoA[$i][$s] eq 'A'){
	$SiteCounts[0] += 2;
      }
      elsif ($InputAoA[$i][$s] eq 'C'){
	$SiteCounts[1] += 2;
      }
      elsif ($InputAoA[$i][$s] eq 'G'){
	$SiteCounts[2] += 2;
      }
      elsif ($InputAoA[$i][$s] eq 'T'){
	$SiteCounts[3] += 2;
      }
      elsif ($InputAoA[$i][$s] eq 'M'){
	$SiteCounts[0]++;
	$SiteCounts[1]++;
      }
      elsif ($InputAoA[$i][$s] eq 'R'){
	$SiteCounts[0]++;
	$SiteCounts[2]++;
      }
      elsif ($InputAoA[$i][$s] eq 'W'){
	$SiteCounts[0]++;
	$SiteCounts[3]++;
      }
      elsif ($InputAoA[$i][$s] eq 'S'){
	$SiteCounts[1]++;
	$SiteCounts[2]++;
      }
      elsif ($InputAoA[$i][$s] eq 'Y'){
	$SiteCounts[1]++;
	$SiteCounts[3]++;
      }
      elsif ($InputAoA[$i][$s] eq 'K'){
	$SiteCounts[2]++;
	$SiteCounts[3]++;
      }      
    }
    print O "$SiteCounts[0]\t$SiteCounts[1]\t$SiteCounts[2]\t$SiteCounts[3]\n";
  }
}

close O;

for ($i = 0; $i < @ChrLines; $i++){
  close $ChrHandles[$i];
}

}


