#!/usr/bin/env perl

use strict;
use warnings;

my $ogsfile = shift;

#get all genes
my %allgenes;
my $fastapath = "/Volumes/LaCie/Projects/Current/insimm/musca/orthopipe/OMA/DB";
opendir(my $fastadir, "$fastapath") || die "Can't open fasta path!\n";
while (readdir $fastadir) {
	open FASTA, "$fastapath/$_" or die "Cannot open $fastapath/$_\n!";
	while (<FASTA>) {
		chomp;
		next unless /^>/;
		my $id = $_;
		$id =~ s/^>(\S+).*?$/$1/;
		$id =~ s/\|/_/g;
		$id =~ s/\s+//g;
		$allgenes{$id}=0
	}
}

#read ogs
open OGS, "$ogsfile" or die;
while (<OGS>) {
	chomp;
	my ($id, $genes) = split(/\t/, $_);
	foreach (split(/,/, $genes)) {
		s/\s+//g;
		s/\|/_/g;
		$allgenes{$_}++;
	}
}

foreach my $gene (keys %allgenes) {
	print "$gene\t$allgenes{$gene}\n" unless $allgenes{$gene} == 1;
}
