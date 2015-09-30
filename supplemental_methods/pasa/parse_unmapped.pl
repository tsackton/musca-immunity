#!/usr/bin/perl

use strict;
use warnings;
use Bio::SeqIO;

my $sp = shift;

#first get list of unmapped reads
my %unmapped;
my $list = "/home/tim/projects/inf_pasa/$sp/compreh_init_build/compreh_init_build.details";
open LST, $list or die;

while (<LST>) {
	chomp;
	my ($gene, $trans, $status) = split;
	next if $status eq "pasa";
	next if $status eq "InvalidQualityAlignment_YES_PASAmap";
	$unmapped{$trans} = $gene;
}

#get length info for transcript

open TRIN, "/home/tim/projects/inf_pasa/$sp/Trinity.$sp.fa.clean" or die;
my %length;

while (<TRIN>) {
	chomp;
	next unless s/^>//;
	#>comp0_c0_seq1 len=340 path=[6:0-263 2810:264-339]
	my ($trans, $len, undef) = split;
	$len =~ s/len=//;
	$length{$trans} = $len;
}

open GFF, "/home/tim/projects/inf_pasa/$sp/Trinity.$sp.fa.clean.transdecoder.gff3" or die;
my %cds;
while (<GFF>) {
	next if /^#/;
	next if /^\s*$/;
	chomp;
	my @fields = split(/\t/, $_);
	next unless $fields[2] eq "CDS";
	my $cdslen = $fields[4]-$fields[3]+1;
	if ($fields[8] =~ /Parent=(\w+)\|/) {
		my $id = $1;
		$cds{$id} = $cdslen;
	}
}

my %keep;

#now read through unmapped list and keep based on heuristics
foreach my $trans (keys %unmapped) {
	my $totlen = $length{$trans};
	my $cdslen = $cds{$trans} || 0;
	my $cdsfrac = $cdslen / $totlen;
	
	if ($cdslen < 180 and $totlen < 200) {
		next;
	}
	elsif ($cdslen < 180 and $totlen > 200) {
		$keep{$trans} = $unmapped{$trans};
	}
	elsif ($cdslen < 450 and $cdsfrac < 0.75) {
		next;
	}
	elsif ($cdslen < 450 and $cdsfrac >= 0.75) {
		$keep{$trans} = $unmapped{$trans};
	}
	elsif ($cdslen >= 450) {
		$keep{$trans} = $unmapped{$trans};
	}
}

#get sequences
my $protseqfile = "/home/tim/projects/inf_pasa/$sp/Trinity.$sp.fa.clean.transdecoder.pep";
my $rnaseqfile = "/home/tim/projects/inf_pasa/$sp/Trinity.$sp.fa.clean";

my $protseq = Bio::SeqIO->new(-file=>$protseqfile);
my $rnaseq = Bio::SeqIO->new(-file=>$rnaseqfile);

my %prots;
my %rnas;

my %skip;

while(my $seq = $protseq->next_seq) {
	my $id = $seq->display_id;
	my $len = $seq->length;
	$id =~ s/\|.*$//;
	my $seqtxt = $seq->seq;
	my $desc = $seq->desc;
	my $type = undef;
	if ($desc =~ /.*type:(\w+).*/) {
		$type = $1;
	}

	if ($type eq "internal" and $len < 200) {
		$skip{$id}++;
		next;
	}

	$prots{$id} = $seqtxt;
}

while(my $seq = $rnaseq->next_seq) {
	my $id = $seq->display_id;
	$id =~ s/\|.*$//;
	my $seqtxt = $seq->seq;
	$rnas{$id} = $seqtxt;
}


open PROT, ">pasa_out/${sp}.unmapped.proteins.fa" or die;
open RNA, ">pasa_out/${sp}.unmapped.transcripts.fa" or die;
open KEY, ">pasa_out/${sp}.unmapped.isokey.txt" or die;

foreach my $trans (keys %keep) {
	next if exists($skip{$trans});
	my $gene = $keep{$trans};
	print KEY "$gene\t$trans\n" if exists($rnas{$trans});
	print PROT ">$trans\n$prots{$trans}\n" if exists($prots{$trans});
	print RNA ">$trans\n$rnas{$trans}\n" if exists($rnas{$trans});
}
	

