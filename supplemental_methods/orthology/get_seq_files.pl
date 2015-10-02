#!/usr/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;

#read ogs
my %ogs;
open OGS, "ogs_updated_hmm.txt" or die;
while (<OGS>) {
	chomp;
	my ($id, $genes) = split(/\t/, $_);
	$ogs{$id} = $genes;
}

#get all genes
my %seqs;
my $fastapath = "/Volumes/LaCie/Projects/Current/insimm/musca/ogs/OMA/DB";
opendir(my $fastadir, "$fastapath") || die "Can't open fasta path!\n";
while (readdir $fastadir) {
	my $allprot = Bio::SeqIO->new(-file=>"$fastapath/$_", -format=>"fasta");
	while (my $seq = $allprot->next_seq) {
		my $id = $seq->display_id;
		my $seqtext = $seq->seq;
		$seqs{$id}=$seqtext;
	}
}
my $counter = 0;
my $ogcount = 0;

foreach my $og (keys %ogs) {
	$ogcount++;
	if ($ogcount > 270) {
		$counter++;
		$ogcount=0;
	}
	my $dir = $counter + 100;
	system("mkdir -p newogs_seq/$dir");
	open OUT, ">./newogs_seq/$dir/$og.fa" or die;
	my @genes = split(/,/, $ogs{$og});
	foreach my $gene (@genes) {
		unless ($seqs{$gene}) {
			warn "Cannot find proteins for $gene\n";
			next;
		}
		print OUT ">$gene\n$seqs{$gene}\n";
	}
}
