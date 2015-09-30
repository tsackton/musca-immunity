#!/usr/bin/env perl

use strict;
use warnings;
use List::Util qw(min max);
use Bio::SeqIO;

#define files
my $outfile = "mdom_strata.txt";
my $reblastfile = "mdom_reblast.txt";
my $nonprotfile = "mdom_nonprot.txt";

#get strata key
open SKEY, "strata_key.txt" or die;
my %spkey;
my %sp_to_strat;
my %strat_key;

while (<SKEY>) {
	chomp;
	next if /^spname/;
	my @fields = split(/\t/, $_);
	$spkey{$fields[1]} = $fields[0];
	my $strattext = $fields[2];
	
	my ($strat_num, $strat_tax) = split(/:/, $strattext);
	$sp_to_strat{$fields[1]} = $strat_num;
	$strat_key{$strat_num} = $strat_tax;
}

my @taxa = map { $strat_key{$_}} sort { $a <=> $b } keys %strat_key;

#get full gene list
open SPKEY, "../annot/mdom.final.protkey" or die;
my %genes;
while (<SPKEY>) {
	chomp;
	s/mdom\|//g;
	my ($gene, $trans) = split(/\t/, $_);
	$genes{$trans} = $gene;
}
my %blast;

#get protein list and species key
my $protfile = Bio::SeqIO->new(-file=>"../annot/mdom.final.proteins.fa", -format=>"fasta");
my %proteins;
while (my $prot = $protfile->next_seq) {
	my $id = $prot->display_id;
	my $len = length($prot->seq);
	$proteins{$id}=$len;
	$blast{$id}{'musDom'} = 0
}

#process blast

my @bfiles = <*.out>;
foreach my $file (@bfiles) {
	my $targetsp = $file;
	$targetsp =~ s/mdom_results_//;
	my $infile;
	if ($file =~ /\.gz/) {
		$infile = "gzip -cd $file | ";	
		$targetsp =~ s/\.out\.gz//;
	}
	else {
		$infile = $file;
		$targetsp =~ s/\.out//;	
	}
	open NM, "$infile" or die;
	while (<NM>) {
		my @fields = split(/\s+/, $_);
		my $q = $fields[0];	
		my $qlen = $proteins{$q};
		my $evalue = $fields[10];
		next unless $fields[3]/$qlen >= 0.40;
		next unless $fields[2] > 20;
		my $species = $targetsp;
		$blast{$q}{$species} = exists($blast{$q}{$species}) ? min($blast{$q}{$species}, $evalue) : $evalue;
	}
}

open OUT, ">$outfile" or die;
print OUT "id\tclass\tstrata\t" . join("\t", @taxa) . "\n";
open RB, ">$reblastfile" or die;
open NC, ">$nonprotfile" or die;

foreach my $id (keys %genes) {
	my $type;
	
	if ($id =~ /^PASA/) {
		$type = "novel";
	}
	elsif ($id =~ /^comp/) {
		$type = "unmapped";
	}
	else {
		$type = "known";
	}
	
	if (!exists($blast{$id})) {
		if (exists($proteins{$id})) {
			print RB "$id\n";
#			$type = "blastfail";
		}
		else {
			print NC "$id\n";
			$type = "noncoding";
		}
		
		my $nohit_strata = "Muscidae";
		my $nohit_fields = "\t0" x scalar(@taxa);
		print OUT "$id\t$type\t$nohit_strata" . "$nohit_fields\n";
		next;
	}
	
	my %spcounts;
	my @species = keys %{$blast{$id}};
	
	foreach (@species) {
		my $stra = $sp_to_strat{$_};
		$spcounts{$stra}++
	}
	
	my $deepest = max(keys %spcounts);
	my $deepest_print = $strat_key{$deepest};
	
	my @taxct;
	foreach my $tx (@taxa) {
		my $count = 0;
		foreach my $spe (@species) {
			$count++ if $strat_key{$sp_to_strat{$spe}} eq $tx;
		}
		push @taxct, $count;	
	}

	print OUT "$id\t$type\t$deepest_print\t" . join("\t", @taxct) . "\n";
	
}
		