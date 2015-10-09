#!/usr/bin/env perl

use strict;
use warnings;

#make OGS file from simulated data, and get dup stats from .rel.txt files

open my $ogsfile, ">sim_ogs.txt" or die;
open my $ogsfilect, ">sim_ogs_ct.txt" or die;

my @ogslabs = qw(musDom  gloMor  droWil  droPse  droYak  droMel  droAna  droVir  droMoj  culQui  aedAeg  anoDar  anoSte  anoGam);
print $ogsfilect "OGS\t" . join("\t", @ogslabs) . "\n";
print $ogsfile "OGS\tgenes\n";

my %dups;

my @rates = qw(0.00057 0.00114 0.002 0.005 0.007 0.009 0.014 0.023 0.034 0.057 0.114 0.341);
foreach my $rate (@rates) {
	my $rateid = $rate;
	$rateid =~ s/^0\.//;
	opendir(my $simdir, $rate) || die;
	while(readdir $simdir) {
		my $file = "$rate/$_";
		next unless $_ =~ /\.pruned\.leafmap/;
		next if $_ =~ /^sim\.1000\./;
		my $simrun = $_;
		$simrun =~ s/\.pruned\.leafmap//;
		my $ogsid = $simrun . "." . $rateid;
		my %ogsct = qw(musDom 0 gloMor 0 droWil 0 droPse 0 droYak 0 droMel 0 droAna 0 droVir 0 droMoj 0 culQui 0 aedAeg 0 anoDar 0 anoSte 0 anoGam 0);
		my @ogsgenes;
		my %spkey;
		open my $infile, $file or die "Cannot open $file\n";
		while (<$infile>) {
			chomp;
			my ($gene, $sp) = split;
			$ogsct{$sp}++;
			push @ogsgenes, "${sp}_${ogsid}_${gene}";
			$spkey{$gene}=$sp;
		}
		print $ogsfile "$ogsid\t" . join(",", @ogsgenes) . "\n";
		print $ogsfilect "$ogsid";
		foreach my $label (@ogslabs) {
			print $ogsfilect "\t$ogsct{$label}";
		}
		print $ogsfilect "\n";
		
		#now read rel file and make dup/loss table 
		open my $rel, "$rate/$simrun.pruned.rel.txt" or die "Cannot open rel file for $rate, $simrun\n";
		my $curogs = $ogsid;
		$dups{$curogs}{'dup'}{'total'}=0;
		$dups{$curogs}{'loss'}{'total'}=0;	
		while (<$rel>) {
			chomp;
			my @fields = split;
			if ($fields[0] eq "gene") {
				my $spec_id = $spkey{$fields[1]};
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
}

open my $dupfile, ">ogs_dup_ct.txt" or die;
my @options = qw(total musDom gloMor 7 droMel droYak 13 12 droAna 11 droPse 10 droWil 9 8 droMoj droVir 6 anoGam anoSte 4 anoDar 3 aedAeg culQui 5 2 1);
print $dupfile "ogs\ttype\t" . join("\t", @options) . "\n";
foreach my $ogs (keys %dups) {
	my $dupline = "$ogs\tdup";
	my $lossline = "$ogs\tloss";
	my $totline = "$ogs\ttotal";
	foreach my $col (@options) {
		my $dup;
		my $loss;
		my $tot;
		$dup = defined($dups{$ogs}{'dup'}->{$col}) ? $dups{$ogs}{'dup'}->{$col} : "NP";
		$loss = defined($dups{$ogs}{'loss'}->{$col}) ? $dups{$ogs}{'loss'}->{$col} : "NP";
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
		$dupline .= "\t$dup";
		$lossline .= "\t$loss";
		$totline .= "\t$tot";
	}
	print $dupfile "$dupline\n$lossline\n$totline\n";
}
	
