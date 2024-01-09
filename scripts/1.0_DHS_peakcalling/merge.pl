#!/bin/perl -w
use strict;
use lib '/home/orchard/Scripts/new_scripts_file/Modules/';
use range;

my $file_of_peaks = shift; ###sorted!!!!!
my $necessary_overlap = shift;

###print the header now
my @header = ("region_id", "region_chromosome", "region_start", "region_end", "peak_chromosome", "peak_start", "peak_end", "srx_id", "peak_name", "peak_score");
print join("\t", @header), "\n";
my $region_id = 0; 

open(PEAKS, $file_of_peaks);

my $current_range;
my @ranges_in_current_range;

while (<PEAKS>) {
	chomp;
	my @line = split /\t/, $_;
	my $range = range->new();
	$range->set_namespace(shift @line);
	$range->set_start(shift @line);
	$range->set_end(shift @line);
	my $range_id = join("\t", @line);
	$range->set_id($range_id); ###"srx_id	peak_name	peak_score"
	#print $range->start(), "\n";
	#print @line;
	#printf STDERR "$range_id\n";

	if (!defined($current_range)) {
		$current_range = $range;
		push @ranges_in_current_range, $range;
	} else {
		my $intersect = 0;
		if ($current_range->overlaps($range)) {
			$intersect = $range->intersect($current_range);
			$intersect = $intersect->get_length();
		}
		
		###if the range overlaps with the "current range", then merge the ranges and continue
		if ($intersect >= $necessary_overlap) {
			$current_range = $current_range->merge($range);
			push @ranges_in_current_range, $range;
		} else {
			print_merged($current_range, @ranges_in_current_range);
			$current_range = $range;
			@ranges_in_current_range = ($range);
		}


		###otherwise process the old range
	}
}

###take care of the last one!
print_merged($current_range, @ranges_in_current_range);

close PEAKS;

sub print_merged {
	my $current_range = shift;
	my @ranges_in_current_range = @_;
	my ($region_chromosome, $region_start, $region_end) = ($current_range->namespace(), $current_range->start(), $current_range->end());
	$region_id++;
	foreach my $range (@ranges_in_current_range) {
		my ($range_chromosome, $range_start, $range_end, $range_id) = ($range->namespace(), $range->start(), $range->end(), $range->id());
		my ($srx_id, $peak_name, $peak_score) = split /\t/, $range_id;
		my @to_print = ($region_id, $region_chromosome, $region_start, $region_end, $range_chromosome, $range_start, $range_end, $srx_id, $peak_name, $peak_score);
		my $to_print = join("\t", @to_print);
		print "$to_print\n";
	}
}
