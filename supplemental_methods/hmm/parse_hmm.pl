#!/usr/bin/env perl

use strict;
use warnings;
use List::Util qw(min);

open OUT, ">hmm_parsed.txt" or die;
open IN, "target_out.tab" or die;
my %targets;
while (<IN>) {
	next if /^#/;
	chomp;
	my @fields = split;
	my $target = $fields[0];
	my $query = $fields[2];
	my $evalue = $fields[4];
	my $onedome = $fields[7];
	
	next if $onedome > 1; #skip cases where the score based on best domain is very weak
	next unless $target =~ /^mDom/; #only get mDom proteins
	$target =~ s/^mDom\|//;
	next if $target =~ /^comp/; #clean up old models
	print OUT "$target\t$query\t$evalue\n";
}