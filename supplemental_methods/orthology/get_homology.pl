#!/usr/local/bin/perl

use strict;
use warnings;

#read through all orth.txt files
#get mdom<->dmel orthologs
#these can be 1:1, 1:many, many:1, or many:many

open MOG, "final_ogs.txt" or die;
my %mogs;
while (<MOG>) {
	chomp;
	my ($ogs, $genes) = split(/\s+/, $_);
	$mogs{$ogs}++;
}


my %homs;
my %type;

my @relfiles = <homology/orth/*.orth.txt>;
foreach my $file (@relfiles) {
	my $ogs = $file;
	$ogs =~ s/homology\/orth\///;
	$ogs =~ s/\.orth\.txt//;
#	print "$ogs\n";
	next unless exists ($mogs{$ogs});
	open REL, $file or die;
	while (<REL>) {
		chomp;
		my @fields = split;
		next unless $fields[0] eq "droMel" and $fields[1] eq "musDom";
		$fields[3] =~ s/musDom_//;
		$fields[2] =~ s/droMel_//;
		push @{$homs{$fields[3]}}, $fields[2];
		$fields[4] = "n" if $fields[4] > 1;
		$fields[5] = "n" if $fields[5] > 1;
		$type{$fields[3]} = "$fields[5]:$fields[4]";
	}
}

open HOM, ">dmel_mdom_homology.txt" or die;
print HOM "mdom\ttype\tdmel\n";

foreach my $mdomid (sort keys %homs) {
	my @dmelids = @{$homs{$mdomid}};
	my $orthtype = $type{$mdomid};
	print HOM "$mdomid\t$orthtype\t" . join(",", @dmelids) . "\n";
}
