#!/usr/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;

my $file = shift;
open my $out, ">$file.lengths" or die;
my $inseq = Bio::SeqIO->new(-file => $file);
while (my $seq=$inseq->next_seq) {
	my $sequence = $seq->seq;
	my $fullheader = $seq->id . " " . $seq->desc;
	my $id = $seq->id;
	my $length = $seq->length;
	print $out "$id\t$length\n";
}
	