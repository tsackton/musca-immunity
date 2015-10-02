#!/usr/bin/env perl

use strict;
use warnings;
use Graph;
use 5.010;

#basic strategy - still a bit of a work in progress
#read in OGS, make graph where each protein is a node
#get the best HMM for each protein (can be none if protein does not hit anything)
#for each protein, add an edge between that protein and all proteins in the best HMM
#make sure to add singleton edges too somehow
#now extract CC from graph, should be exactly right

#read OG file and get the OG that each gene belongs to (or none)
open OGS, "ogs.txt" or die;
my %gene_to_ogs;
my %ogs_to_members;
my $newogs = Graph::Undirected->new;

while (<OGS>) {
	chomp;
	my ($ogs, $genes) = split(/\t/, $_);
	$genes =~ s/^\s+//;
	$genes =~ s/\s+$//;
	my @genes = split(/,/, $genes);
	foreach my $gene (@genes) {
		$gene_to_ogs{$gene} = $ogs;
		foreach my $gene2 (@genes) {
			$newogs->add_edge($gene, $gene2); #treats orthogroups as connected graphs
		}
	}
	$ogs_to_members{$ogs} = \@genes;
}

#get full list of proteins and assign to "none" any without an OGS assignment
my $fastapath = "/Volumes/LaCie/Projects/Current/insimm/musca/orthopipe/OMA/DB";
opendir(my $fastadir, "$fastapath") || die "Can't open fasta path!\n";
while (readdir $fastadir) {
	open FASTA, "$fastapath/$_" or die "Cannot open $fastapath/$_\n!";
	while (<FASTA>) {
		chomp;
		next unless /^>/;
		my $id = $_;
		$id =~ s/^>(\S+).*?$/$1/;
		$gene_to_ogs{$id} = "none" unless exists ($gene_to_ogs{$id});
		$newogs->add_edge($id, $id) if $gene_to_ogs{$id} eq "none";
	}
}

#read HMM file and get best OG for each gene (or none)
my $hmmcmd = "gzip -cd /Volumes/LaCie/Projects/Current/insimm/musca/orthopipe/hmmscan/hits/merged-hits.txt.gz | ";
open HMM, "$hmmcmd" or die;
my %besthmm;
while (<HMM>) {
	chomp;
	next if /^#/;
	my @fields = split;
	my $prot = $fields[0];
	my $hmm = $fields[2];
	my $evalue = $fields[4];
	if (exists $besthmm{$prot}{'hmm'}) {
		#debugging -- see what is actually being compared
		print "$prot: testing if $besthmm{$prot}{'evalue'} > $evalue\n";
		if ($besthmm{$prot}{'evalue'} > $evalue) { #this is probably not ideal as it means that a gene with 2 or more 0s will get assigned a best hmm based on the order the hits appear in the list
			#fix this!! should have everything that is equal as one set
			print "\t$evalue is better than $besthmm{$prot}{'evalue'}\n";
			$besthmm{$prot}{'hmm'}=$hmm;
			$besthmm{$prot}{'evalue'}=$evalue;
			next;
		}
		else {
			next;
		}
	}

	$besthmm{$prot}{'hmm'}=$hmm;
	$besthmm{$prot}{'evalue'}=$evalue;
}

#now merge duplicates based on hmmsearch (best HMM for each protein)
foreach my $protid (keys %besthmm) {
	my $besthit = $besthmm{$protid}{'hmm'};
	my $oldogs = $gene_to_ogs{$protid} || "NA";
	print "$protid\t$oldogs\t$besthit\n";
	if ($besthit eq $oldogs) {
		next;
	}
	else {
		my @members = @{$ogs_to_members{$besthit}};
		foreach my $member (@members) {
			$newogs->add_edge($member, $protid);
		}
	}
}

#extract connected components

open NEW, ">ogs_updated_hmm_2.txt" or die;
my @cc = $newogs->connected_components();

my $newid = 0;
foreach my $component (@cc) {
#	print "$component\n";
	#each component gets its own newid
	$newid++;
	my @ids = @{$component};
#	print "@ids\n";
	print NEW "$newid\t" . join(",", @ids) . "\n";
}
