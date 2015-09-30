#!/usr/bin/perl

use strict;
use warnings;

my $oldgff = shift;
my $newgff = "$oldgff.updated";

open IN, "$oldgff" or die;
open OUT, ">$newgff" or die;

my %seen_ids;
my %subs;

while (<IN>) {
	if (/^#/) {
		print OUT;
		next;
	}

	if (/^\s*$/) {
		print OUT;
		next;
	}

	chomp;
	my @fields = split(/\t/, $_);
	my %notes = split(/[=;]/, $fields[8]);

	#get id and parent
	my $id = $notes{'ID'};
	my $parent = $notes{'Parent'} || "NONE";

	#if ID already exists, update it and indicate that any further parents with that id will be replaced
	if (exists($seen_ids{$id})) {
		my $counter = $seen_ids{$id}+1;
		my $origid = $id;
		$id = "$id.$counter";
		$subs{$origid} = $id;
	}

	if (exists($subs{$parent})) {
		$parent = $subs{$parent};
	}

	$seen_ids{$id}++;

	if ($parent ne "NONE" ) {
		$fields[8] = "ID=$id;Parent=$parent;";
	}
	else {
		$fields[8] = "ID=$id;";
	}

	my $printline = join("\t", @fields);
	print OUT "$printline\n";

}
