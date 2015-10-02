#!/usr/bin/env perl

use strict;
use warnings;
use 5.012;

my %ogs;

#make OGS file from HOGFasta input
my $OMApath = "/Volumes/LaCie/Projects/Current/insimm/musca/ogs/OMA";

opendir(my $hogpath, "$OMApath/Output/HOGFasta") || die "Can't open HOGFasta path!\n";
while (readdir $hogpath) {
	my $og = $_;
	$og =~ s/\.fa//;
	open FASTA, "$OMApath/Output/HOGFasta/$_" or die "Cannot open $hogpath/$_\n!";
	while (<FASTA>) {
		chomp;
		next unless /^>/;
		my $id = $_;
		$id =~ s/>(\S+)\s+.*$/$1/;
		push @{$ogs{$og}}, $id;
	}
}

open OUT, ">ogs.txt" or die;

foreach my $og (keys %ogs) {
	print OUT "$og\t" . join(",", @{$ogs{$og}}) . "\n";
}
		