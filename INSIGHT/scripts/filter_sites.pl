#!/usr/bin/perl -w

$in1 = $ARGV[0];
$in2 = $ARGV[1];
$out = $ARGV[2];
print `date`;

chomp($in1);
chomp($in2);
%hash=();

open(TEXT,$in1) or die "Cant open file\n";
while (<TEXT>)
{
	chomp($_);
	@arr = split("\t", $_);
	for($i=$arr[1];$i<=$arr[2];$i++){
		$hash{$arr[0]}{$i} = 0; 
	}
}
close TEXT;
print `date`;

open(IN,$in2) or die "Cant open file\n";
open(OUT,">$out") or die "Cant open file\n";
while (<IN>)
{
	chomp($_);
	if($_ =~ m/^\#/){
		$flag=0;
		$l = $_;
	}
	else{
		@arr = split(/\s+/, $_);
		
		if(exists $hash{$arr[0]}{$arr[1]}){
		if($flag==0){print OUT $l,"\n";$flag=1;}
			print OUT $_,"\n";
		}
	}
}
close IN;
close OUT;
print `date`;
