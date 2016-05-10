#!/usr/bin/env perl

use strict;
use warnings;
use List::Util qw(min);

my $file = shift;
my $species = $file;
print STDERR "$species\n";
$species =~ s/^raw\/(.+?)\.fasta\.tab$/$1/;
print STDERR "$species\n";

open IN, "$file" or die;
my %targets;
while (<IN>) {
	next if /^#/;
	chomp;
	my @fields = split;
	my $target = $fields[0];
	my $query = $fields[2];
	my $evalue = $fields[4];
	my $onedome = $fields[7];
	
	next if $onedome > 0.001; #skip cases where the score based on best domain is very weak
	next if $evalue > 1e-5; #skip cases without strong evidence 
	
	if ($species eq "mdom.final.proteins") {
		chop($target);
	}
	
	print "$species\t$target\t$query\t$evalue\n";
}