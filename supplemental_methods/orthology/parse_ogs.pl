#!/usr/bin/perl

use strict;
use warnings;

#make two keys: ogs->species counts and mdom->ogs

open IN, "final_ogs.txt" or die;
my %mdom;
my %dmel;
my @spkey = qw(musDom gloMor droMel droYak droAna droPse droWil droMoj droVir anoGam anoDar anoSte aedAeg culQui);

open OGS, ">ogs_stats.txt" or die;
print OGS "ogs\t" . join("\t", @spkey) . "\n";
while (<IN>) {
	chomp;
	my ($ogs, $genes) = split(/\s+/, $_);
	my @genes = split(/,/, $genes);
	my %sp;
	foreach my $gene (@genes) {
		my ($sp, $geneid) = split(/[_\|]/, $gene, 2);
		$sp{$sp}++;
		die "$geneid in multiple OGS rows!" if exists($mdom{$geneid});
		$mdom{$geneid} = $ogs if $sp eq "musDom";
		$dmel{$geneid} = $ogs if $sp eq "droMel";
	}
	print OGS "$ogs"; 
	foreach my $spec (@spkey) {
		my $num = $sp{$spec} || 0;
		print OGS "\t$num";
	}
	print OGS "\n";
}

open MDOM, ">mdom_to_ogs.txt" or die;
print MDOM "mdomid\togsid\n";

open DMEL, ">dmel_to_ogs.txt" or die;
print DMEL "dmelid\togsid\n";

foreach my $mdomid (keys %mdom) {
	print MDOM "$mdomid\t$mdom{$mdomid}\n";
}

foreach my $dmelid (keys %dmel) {
	print DMEL "$dmelid\t$dmel{$dmelid}\n";
}