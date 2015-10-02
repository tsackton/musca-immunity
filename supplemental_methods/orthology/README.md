Orthology pipeline
==================


Update based on HMM search
-----------

1. Run OMA (default options)

2. Make ogs.txt file from OMA output (make_ogs_from_HOG.pl)

3. Make alignments from HOGs (build_hmms_from_HOGS.sh)

4. Search hmms (search_hmms.sh)

5. Update OGS after hmmscan (fix_ogs_after_hmmscan.pl)

6. Make new sequence files (get_seq_files.pl)

7. Check OGS after HMM update (check_ogs.pl)


Tree-building
--------

1. Align new OGS files (align_ogs.sh)

2. Run raxml (run_trees.sh)


Update OGS based on trees
------

1. Run treesplit_init.pl (changing input / output dirs as needed)

2. Run treefix (run_treefix.sh)

3. Rerun tree splitting (treesplit_final.pl)

4. Rerun treefix on changed OGS after final treesplit (run_treefix.sh)


Final processessing
--------

1. Check ogs again (check_ogs.pl)

2. Annotate trees (annotate_trees.sh)

3. Parse and get stats (parse_ogs.pl, get_dupstats.pl, get_homology.pl)
