#!/bin/perl -w
use strict;
use Getopt::Long;

=head1 NAME

filter.pl

=head1 SYNOPSIS

=head1 DESCRIPTION

A general utility used for quickly filtering files that I wrote while working on my Master's thesis.

=head1 AUTHOR

=head1 APPENDIX

What follows is a description of the modes that may be passed as the '-m' parameter.

=head2 Mode 1

 Required arguments:
	-m
	-i
	-l
	-f

	Given a file, a list of items, and a column number (or list of column numbers), checks the columns in each line against the items in the file.  It the columns in the row match one of the items in the file, the row is printed; otherwise, it is not printed.  the -v=1 flag inverts the print/don't print action. Columns are indexed starting from 1. The list of items should use the same delimiter as the file to be filtered.

	Obama	United_states	president
	Merkel	Germany			chancellor
	Biden	United_states	vice_president
	Tusk	Poland			prime_minister
	Brown	England			prime_minister
	
	if -l passes the following file:

	United_states
	Poland

	and -f is set to 2 (second columm), then the following is printed:

	Obama	United_states	president
	Biden	United_states	vice_president
	Tusk	Poland			prime_minister

	if -v were passed as a flag, the inverse would be printed


=head2 Mode 2

	Accepts a list of "items to keep", a list of "items to exclude", and a file to filter.  Scans through the file to filter, and uses the "items to keep" or "items not to keep" as sort of "print/don't print" switches.  For each line in the file to filter, the script checks whether or not the line is contained in one of the lists.  If it is in the "items to keep list", the switch is set to "print", and the line is printed.  If the item is in the "not to keep" list, the switch is set to "don't print", and the line is not printed.  If a line does not appear in either of the lists, it is printed/not printed based on the state of the switch.  The switch starts in the "print" position.


=cut

###how do you pass a tab as a command line argument??

if(scalar @ARGV == 0){
	my $usage = <<'END_USAGE';

	Usage:
		-m mode (required)
		-i input_file
		-d delimiter (optional; defaults to tab)
		-l list_of_items (to be used with mode 1)
		-f fields_by_which_to_filter (to be used with mode 1; one or more integer values, comma-separated)
		-le list_of_items_to_exclude (to be used with mode 2)
		-lk list_of_items_to_keep (to be used with mode 2)
		-v (flag for inverse; optional)
		-h header_is_present (flag indicating the presence of a header line in input_file; prints the first line of input_file without filtering)

	Modes:
		1: Filter according to one or more columns
		2: Filter according to entire lines

END_USAGE

	print $usage;
	exit;
}

###get options
my ($mode, $delimiter, $list_of_items, $list_of_items_to_keep, $list_of_items_to_exclude, $input_file, $inverse, $header, @fields);
GetOptions('m=i' => \$mode, 'i=s' => \$input_file, 'd:s' => \$delimiter, 'l:s' => \$list_of_items, 'le:s' => \$list_of_items_to_exclude, 'lk:s' => \$list_of_items_to_keep, 'v' => \$inverse, 'f:s' => \@fields, 'h' => \$header);
$mode || die "No mode option provided\n";
$input_file || die "No input file provided\n";
$delimiter || ($delimiter = "\t");
@fields = split(/,/, join(',',@fields)) if @fields;


if ($mode == 1){
	$list_of_items || die "No list of items provided\n";
	$inverse || ($inverse = 0);
	mode_1($inverse, \@fields, $delimiter, $list_of_items);
} elsif ($mode == 2){
	$list_of_items_to_exclude || die "No list of items to exclude\n";
	$list_of_items_to_keep || die "No list of items to keep\n";
	mode_2($list_of_items_to_keep, $list_of_items_to_exclude);
}


################# SUBFUNCTIONS BELOW THIS LINE ###################

sub mode_1 {
	my ($inverse, $fields, $delimiter, $list_of_items) = @_;
	my @fields = @{$fields};
	my %to_keep;

	open(KEEP, $list_of_items);

	while(<KEEP>){
		chomp;
		$to_keep{$_} = 1;
	}

	close KEEP;
	
	open(FILTER, $input_file);
	
	while(<FILTER>){
		chomp;
		if ($. == 1 && $header) {
			print "$_\n";
			next;
		}
		my @line = split /$delimiter/, $_;
		my @to_check;
		foreach my $field (@fields){
			push @to_check, $line[$field-1];
		}
		my $to_check = join "$delimiter", @to_check;
		if (exists $to_keep{$to_check}){
			print "$_\n" if $inverse == 0;
		} else {
			print "$_\n" if $inverse == 1;
		}
	}
	
}

sub mode_2 {
	my ($list_of_items_to_keep, $list_of_items_to_exclude) = @_;
	my (%items_to_keep, %items_to_exclude);
	my $print_switch = "yes";
	open(KEEP, $list_of_items_to_keep);

	while (<KEEP>){
		chomp;
		if ($. == 1 && $header){
			print "$_\n";
			next;
		}
		$items_to_keep{$_} = 1;
	}

	close KEEP;

	open(EXCLUDE, $list_of_items_to_exclude);

	while (<EXCLUDE>){
		chomp;
		$items_to_exclude{$_} = 1;
	}

	close EXCLUDE;

	open(FILTER, $input_file);

	while(<FILTER>){
		chomp;
		if ( exists($items_to_keep{$_}) ) {
			$print_switch = "yes";
		} elsif ( exists($items_to_exclude{$_}) ) {
			$print_switch = "no";
		}

		if ($print_switch eq "yes"){
			print "$_\n";
		}
	}

	close FILTER;
}
