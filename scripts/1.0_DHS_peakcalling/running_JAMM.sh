#!/bin/bash

directory_with_sample_files=$1
output_directory=$2  ##this needs to be premade
resolution=$3 #peak, region, window
number_processors_to_use=$4

genome_size_file="/home/orchard/General/hg19_chromosome_size_file.txt"

bash /home/orchard/jamm/JAMM1.0.7rev1/JAMM.sh -g $genome_size_file -f 1 -d y -s $directory_with_sample_files -r $resolution -o $output_directory -p $number_processors_to_use

#OPTIONS for JAMM:
#   -s      directory containing Sample files (required) #argument 1
#   -g      Genome size file (required) #already set
#   -o      Output directory (required) #argument 2
#   -c      directory containing input or Control files #don't have
#   -f      Fragment length(s) (default: estimated) # Will use -f 1, because this is mentioned in the paper (table 1)
#   -r      Resolution, peak or region or window (default: peak) #argument 3
#   -m      Mode, normal or narrow (default: normal) #will leave at the default
#   -i      clustering Initialization window selection, deterministic or stochastic (default: deterministic)
#   -b	   Bin Size (default: estimated)
#   -w      minimum Window size (default: 2 --- Note: this means minimum_window_size = bin_size x the_value_of_-w)
#   -e	   window Enrichment cutoff, auto or any numeric value (default: 1 --- Set this to "auto" to estimate the window enrichment cutoff)
#   -d	   keep PCR Dupicates in single-end mode, y or n (default: n --- if -t is "paired", this option has no effect) #we've already removed them I think
#   -t	   Type, single or paired (default: single, requires BED files. paired requires BEDPE files)
#   -p	   Number of processors used by R scripts (default: 1)


