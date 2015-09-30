#!/usr/bin/perl

my $in = shift;
my $out = "$in.new";

open IN, "$in" or die;
open OUT, ">$out" or die;

while (<IN>) {
	my @fields = split(/\t/, $_);
	s/ID=.*?;// if $fields[2] eq "CDS" or $fields[2] eq "exon";
	print OUT;
}
