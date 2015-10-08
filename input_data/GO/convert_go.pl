#!/usr/bin/perl

use strict;
use warnings;

open IN, "annot_Seqs_20140317_1550.txt" or die;
my %terms;
while (<IN>) {
	chomp;
	if (/^ref\|(XP_\d+)\|.*(GO:\d+)/) {
		push @{$terms{$1}}, $2;
	}
}

open OUT, ">mdom_go_list" or die;

foreach my $xp (keys %terms) {
	foreach my $term (@{$terms{$xp}}) {
		print OUT "$term\tIEA\t$xp\n";
	}
}
