#!/bin/perl -w
use strict;

###after running peakMerge.pl, can use mergedLabels.pl to classify each region according to tissue specificity

#################################################################################################################################################################
#																																								#
#																******SCRIPT INFORMATION******																	#
# Author:																																						#
# Peter Orchard																																					#
# 																																								#
# Last edit:																																					#
# 																																								#
#																																								#
# Usage: 																																						#
# perl peakMerge.pl sorted.DNaseCombine.output.bed toInclude.txt																								#
#																																								#
# Purpose:																																						#
# Merges similar peaks from different DNase-seq experiments, as comparable peaks from these experiments often will differ in their exact start/end point.		#
# If two peaks overlap at all (i.e., the end base of one peak is the start base for another peak), these peaks are considered to be in the same					#
# "peak region". For the SRX experiments of interest, a "peak_regions_file" will be created. For each peak from the relevant experiments, a region will 		#
# be given (chromosome, region_start point, and region_end point). The peak_score (e.g., "36") and peak_range (e.g., "99150-99750" if the peak stretches		#
# from 99150 to 99750 on chromosome 3) are given, as is the peak_coverage (the width of the peak divided by the width of the peak_region--e.g., 80.5 for 		#
# 80.5%).																																						#
#																																								#
# Input files:																																					#
# 1. DNaseCombine.pl output. Should be formatted as follows:																									#
#		chr		peak_start	peak_end	SRX_pk#		max score	strand (always +)																				#
#	 input file must also be sorted by chromosome and start/end (sort -k1,1 -k2n,2 filename.bed > sorted.filename.bed)											#
# 2. "to include" file.  Simple list of SRX experiments to be included in the merging. Should be just one experiment per line.									#
#																																								#
#																																								#
# Output:																																						#
# Outputs the contents of a peak_regions file.  For each peak in the input file (for the experiments that were included), one line should be generated. Format:	#
# chromosome	region_start	region_end	peak	peak_score	peak_range	peak_coverage																		#
#																																								#
#################################################################################################################################################################	


my $bed_file = shift @ARGV;
my $to_include = shift @ARGV;
my %include;
if($to_include ne "all"){
  open(INCL, "$to_include");
  while(<INCL>){
    chomp;
    $include{$_} = 1;
  }
}


my $current_number_of_peaks;
my $current_region_start = 0;
my $current_region_end = 1;

my %current_region;

my $old_chrom; #holds the chromosome label for the current region
my $new_chrom; #holds the chromosome label for the chromosome on the current line

open(BED, $bed_file);

LINE: while(<BED>){
	chomp;
	my @line = split/\t/, $_;
	my $peak_name = $line[3];
	my $tissue = $line[10];
	$peak_name = $peak_name . "::" . $tissue;
	if ($to_include ne "all"){
	  my $peak_name2 = $peak_name;
	  $peak_name2 =~ s/_pk.*//;
	  next LINE if (!exists $include{$peak_name2});
	}
	$new_chrom = $line[0];
	my $peak_start = $line[1];
	my $peak_end = $line[2];
	my $peak_strength = $line[6];
	my $peak_strand = $line[5];
	#check to see if the peak is overlapping with past peaks:
	if (peakStart_in_current_peak_region($peak_start, $current_region_start, $current_region_end, $new_chrom, $old_chrom)){ ##this new peak overlaps with the current peaked region
		if ($peak_end > $current_region_end) { #extend the overall region if the specific peak overlaps with the region end
			$current_region_end = $peak_end;
		}
		##record the new peak
		$current_region{$peak_name} = [$peak_strength, $peak_start, $peak_end];
	} else { #entering a 0 peak region. need to record all the peaks thus far, and reset the "counters"
		my $region_length = $current_region_end - $current_region_start;
		foreach my $peak (sort keys %current_region) {
			my @peakinfo = @{$current_region{$peak}};
			my $strength = $peakinfo[0];
			my $range = $peakinfo[1] . "-" . $peakinfo[2];
			my $peak_length = $peakinfo[2] - $peakinfo[1];
			my $percent_coverage = ($peak_length / $region_length) * 100;
			$percent_coverage = sprintf "%.1f", $percent_coverage;
			print "$old_chrom\t$current_region_start\t$current_region_end\t$peak\t$strength\t$range\t$percent_coverage\n";
		}
		$current_region_start = $peak_start;
		$current_region_end = $peak_end;
		$old_chrom = $new_chrom;
		%current_region = ();
		$current_region{$peak_name} = [$peak_strength, $peak_start, $peak_end];
	}
}

#need to take care of the last region, which has been missed (double check this??):

my $region_length = $current_region_end - $current_region_start;
foreach my $peak (sort keys %current_region) {
	my @peakinfo = @{$current_region{$peak}};
	my $strength = $peakinfo[0];
	my $range = $peakinfo[1] . "-" . $peakinfo[2];
	my $peak_length = $peakinfo[2] - $peakinfo[1];
	my $percent_coverage = ($peak_length / $region_length) * 100;
	$percent_coverage = sprintf "%.1f", $percent_coverage;
	print "$old_chrom\t$current_region_start\t$current_region_end\t$peak\t$strength\t$range\t$percent_coverage\n";
}


sub peakStart_in_current_peak_region {
	my ($peakstart, $current_peak_start, $current_peak_end, $newChrom, $oldChrom) = @_;
	if (($peakstart >= $current_peak_start) && ($peakstart < $current_peak_end) && ($newChrom eq $oldChrom)) {
		return 1;
	} else {
		return 0;
	}
}
