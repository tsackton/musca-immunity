#!/usr/local/bin/perl

use strict;
use warnings;

open MOG, "final_ogs.txt" or die;
my %mogs;
while (<MOG>) {
	chomp;
	my ($ogs, $genes) = split(/\s+/, $_);
	$mogs{$ogs} = $genes;
}

#read through all rel.txt files
#count total dup+losses, mdom-branch dups, dros-branch dups, mosq dups

my %dups;

my @relfiles = <duplication/*.rel.txt>;
foreach my $file (@relfiles) {
	open REL, $file or die;
	my $curogs = $file;
	$curogs =~ s/duplication\///;
	$curogs =~ s/\.rel\.txt//;
	next unless exists($mogs{$curogs});
	$dups{$curogs}{'dup'}{'total'}=0;
	$dups{$curogs}{'loss'}{'total'}=0;	
	while (<REL>) {
		chomp;
		s/\|/_/g;
		my @fields = split;
		if ($fields[0] eq "gene") {
			my ($spec_id, undef) = split(/_/, $fields[1]);
			$dups{$curogs}{'dup'}{$spec_id} = $dups{$curogs}{'dup'}{$spec_id} || 0;
			$dups{$curogs}{'loss'}{$spec_id} = $dups{$curogs}{'loss'}{$spec_id} || 0;
		}
		elsif ($fields[0] eq "spec") {
			$dups{$curogs}{'dup'}{$fields[3]} = $dups{$curogs}{'dup'}{$fields[3]} || 0;
			$dups{$curogs}{'loss'}{$fields[3]} = $dups{$curogs}{'loss'}{$fields[3]} || 0;
		}			
		elsif ($fields[0] eq "dup") {
			unless ($fields[3]) {
				$fields[3] = $fields[2];
			}
			$dups{$curogs}{'dup'}{$fields[3]}++;
			$dups{$curogs}{'dup'}{'total'}++;
		}
		elsif ($fields[0] eq "loss") {
			$dups{$curogs}{'loss'}{$fields[2]}++;
			$dups{$curogs}{'loss'}{'total'}++;
		}
	}
}

open DUPS, ">ogs_dup_ct.txt" or die;
my @options = qw(total musDom gloMor 7 droMel droYak 13 12 droAna 11 droPse 10 droWil 9 8 droMoj droVir 6 anoGam anoSte 4 anoDar 3 aedAeg culQui 5 2 1);
print DUPS "ogs\ttype\t" . join("\t", @options) . "\n";
foreach my $mog (keys %mogs) {
	my $dupline = "$mog\tdup";
	my $lossline = "$mog\tloss";
	my $totline = "$mog\ttotal";
	foreach my $col (@options) {
		my $dup;
		my $loss;
		my $tot;
		if (!exists($dups{$mog})) {
			if ($col =~ m/^[A-Za-z]+$/ and $col =~ m/$mogs{$mog}/) {
				$dup = 0;
				$loss = 0;
				$tot = 0;
			}
			else {
				$dup = "NP";
				$loss = "NP";
				$tot = "NP";
			}
		}
		else {
			$dup = defined($dups{$mog}{'dup'}->{$col}) ? $dups{$mog}{'dup'}->{$col} : "NP";
			$loss = defined($dups{$mog}{'loss'}->{$col}) ? $dups{$mog}{'loss'}->{$col} : "NP";
			if ($dup eq "NP" and $loss eq "NP") {
				$tot = "NP";
			}
			elsif ($dup eq "NP" and $loss ne "NP") {
				$tot = $loss;
			}
			elsif ($loss eq "NP" and $dup ne "NP") {
				$tot = $dup;
			}
			else {
				$tot = $dup + $loss;
			}
		}
		$dupline .= "\t$dup";
		$lossline .= "\t$loss";
		$totline .= "\t$tot";
	}
	print DUPS "$dupline\n$lossline\n$totline\n";
}
	
