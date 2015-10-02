#!/usr/bin/env perl

use strict;
use warnings;

my %models;
my @files = <raxml/*/RAxML_info.*>;

foreach my $file (@files) {
	open IN, $file or die;
	my $model = "NONE";
	while (<IN>) {
		chomp;
		if (/best-scoring AA model.*\s+(\w+)\s+likelihood/) {
			$model = $1;
		}
	}
	$models{$model}++;
}

foreach my $model (keys %models) {
	print "$model\t$models{$model}\n";
}