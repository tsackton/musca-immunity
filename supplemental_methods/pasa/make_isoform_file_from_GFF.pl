#!/usr/bin/perl

use strict;
use warnings;

while (<>) {
	next if /^#/;
	next if /^\s*$/;
	chomp;
	my @fields = split(/\t/, $_);
	next unless $fields[2] eq "mRNA";
	my %notes = split(/[=;]/, $fields[8]);
	next unless (exists($notes{'ID'}) and exists($notes{'Parent'}));
	print "$notes{'Parent'}\t$notes{'ID'}\n";
}
