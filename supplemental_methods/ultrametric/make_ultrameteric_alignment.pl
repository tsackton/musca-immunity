#!/usr/local/bin/perl

use strict;
use warnings;
use Bio::SeqIO;

open LIST, "ogs_for_ultrametric.txt";
my %files;
while (<LIST>) {
	chomp;
	$files{$_}++;
}

my %merged;

foreach my $filename (keys %files) {
	my $seqfile = "../treefix_final/aligns/$filename.in.fa";
	unless (-r $seqfile) {
		warn "Cannot find alignment for $filename\n";
		next;
	}
	my $seqin = Bio::SeqIO->new(-file=>$seqfile, -format=>'fasta');
	while (my $seq = $seqin->next_seq) {
		my $id = $seq->display_id;
		my $seqtext = $seq->seq;
		$id =~ s/([A-Za-z]+)_.*/$1/;
		$merged{$id} .= $seqtext;
	}
}

open OUT, ">ultrametric.fa" or die;

foreach my $sp (keys %merged) {
	print OUT ">$sp\n$merged{$sp}\n";
}
		
		