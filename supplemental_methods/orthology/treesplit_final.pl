#!/usr/bin/env perl

use strict;
use warnings;
use Bio::TreeIO;
use File::Slurp;
use IO::String;

my $treeDir = $ARGV[0];
my @infiles = <$treeDir/*.out.nwk>;
	
open OUT, ">final_ogs_posttreefix.txt" or die;
my $tree_out_dir = "./treeSplit/";
system("mkdir -p $tree_out_dir");

foreach my $infile (@infiles) {
	print "Starting $infile...";
	#if $infile does not end with ;, add one
	#also set flag for low species number
	my $lowsp =0;
	my $treestring = read_file($infile);
#	print "*******Got tree*********\n$treestring\n";
	chomp($treestring);
	$treestring=~s/\n//g;
	$treestring=~s/\s+//g;
	my $spcount = () = $treestring =~ /[_\|]/g;
	$spcount = $spcount/2;
#	print "$treestring\n$spcount\n";
	$lowsp++ if ($spcount < 4);
	my $singleton = 0;
	$singleton++ if ($spcount == 1);
	warn "Got tree with no taxa at $infile!\n" if $spcount < 1;
#	open FIXTREE, ">$infile" or die;
#	if (substr($treestring, -1) eq ";") {
#		print FIXTREE "$treestring\n";
#	}
#	else {
#		print FIXTREE "$treestring;\n";
#	}
	
	my $ogs = $infile;
	$ogs =~ s/^new\/(\d+.*)\.out\.nwk$/$1/;
	print "OGS is $ogs\n";
	
	if ($singleton) {
		system("cp $infile $tree_out_dir/$ogs.final.nwk");
		(my $sing_id = $treestring) =~ s/[\(\);]//g;
		print OUT "$ogs\t$sing_id\n";
		next;
	}
	elsif ($lowsp) {
		system("cp $infile $tree_out_dir/$ogs.final.nwk");
		(my $ogslist = $treestring) =~ s/[\(\);]//g;
		print OUT "$ogs\t$ogslist\n";
		next;
	}
	
	my @trees_to_parse;
	my @final_trees;

	my %sp_to_clade = qw(droMel Fly droAna Fly droYak Fly droPse Fly droMoj Fly droWil Fly droVir Fly gloMor Fly musDom Fly anoSte Mos anoGam Mos anoDar Mos culQui Mos aedAeg Mos);

	#ugly hack to avoid having to write re-rooting code in perl, and to make Tree object more sensible (since it is by def rooted)
#	system("nw_reroot $infile > temp.tree") == 0 or warn "Reroot failed at $infile!\n";
#	my $rooted = "temp.tree";
#	unless (-s $rooted) {
#		warn "Attempting to open a null tree at $infile!\n";
#		system("mv $infile needBranchLengths");
#		next;
#	}		

	#get tree string
	my $treein = read_file($infile);
	my $io = IO::String->new($treein);
#	print "********Got tree**********\n$treein\n";
	my $treeio = Bio::TreeIO->new(-format => "newick", -fh => $io);
	while (my $origtree = $treeio->next_tree ) {
		$origtree->id($ogs);
		$origtree->nodelete(1);
		push @trees_to_parse, $origtree;
	}

	#now check each tree in @trees_to_parse
	while (@trees_to_parse) {
		my $current_tree = shift @trees_to_parse; #get tree to look at
		#now need to make subtrees
		#get root node
		my $root = $current_tree->get_root_node();
		die "Not rooted" unless defined($root);
#		print "$root\n";
#		print "Processing ", $current_tree->id(), "\n";
		#get descendants of root nodes
		my @subnodes = $root->each_Descendent();

		die "Not bifurcating!" unless scalar @subnodes == 2;
		
		my %check_nodes;
		my $round = 1;
		foreach my $node_to_test (@subnodes) {
#			print "Processing node ", $round, "\n";
			my %ids;
			my @nodes_to_check = $node_to_test->get_all_Descendents();
			foreach my $check (@nodes_to_check) {
				if ($check->is_Leaf()) {
					my ($sp, $gene, $prot) = split(/[_\|]/, $check->id(),3);
					my $species = $sp_to_clade{$sp};
#					print "$species\n";
					$ids{$species}++;
				}
			}
			$check_nodes{$node_to_test} = scalar (keys %ids);
			$round++;
		}
	
		#if all check node keys equal 2, then make subtrees and add to @trees_to_parse
		#otherwise, add current_tree to @final_trees
		my $num_two = 0;
		my $num_keys = scalar keys %check_nodes;
		foreach my $key (keys %check_nodes) {
			$num_two++ if $check_nodes{$key} == 2;
			warn "Error with sp/clade assignment found at $infile!\n" if $check_nodes{$key} > 2;
		}
		
#		print "$infile: Number of 2s: $num_two\t Number of Keys: $num_keys\n";
		if ($num_two == $num_keys) {
			#make substree and add to @trees_to_parse
			my $iter = 1;
			foreach my $new_root (@subnodes) {
#				print "$new_root\n";
#				print "$current_tree\n";
#				print "$infile\n";
				my @testnodes = $new_root->each_Descendent();
#				print "@testnodes", "\n";
				my $newtree = Bio::Tree::Tree->new(-root=>$new_root, -nodelete=>1);
				my $newtree_root = $newtree->get_root_node();
				my @testnodes2 = $newtree_root->each_Descendent();
#				print "$newtree\n";
#				print "$newtree_root\n";
#				print "@testnodes2", "\n";
				my $oldid = $current_tree->id();
				my $newid = "$oldid.$iter";
				$newtree->id($newid);
				push @trees_to_parse, $newtree;
#				my $testout = Bio::TreeIO->new(-format=>'newick', -file=>">$newid.nwk");
#				$testout->write_tree($newtree);
				$iter++;
			}
		}
		else {
			push @final_trees, $current_tree;
		}
	}

	foreach my $tree (@final_trees) {
		my $id = $tree->id();
		my $root = $tree->get_root_node();
		my @nodes = $root->get_all_Descendents();
		my @members;
		foreach my $node (@nodes) {
			if ($node->is_Leaf()) {
				push @members, $node->id();
			}
		}
		print OUT "$id\t", join(",", @members), "\n";
		my $final_tree_file = "$id.final.nwk";
		my $treeout = Bio::TreeIO->new(-format => 'newick', -file => ">$tree_out_dir$final_tree_file");
		$treeout->write_tree($tree);
	}
}

#system("cat final_ogs_multi.txt final_ogs_low.txt > final_ogs.txt");