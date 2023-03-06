#!/bin/bash

clang -O3 -o test_xseed.out test_xseed.c xdrop_aligner.c ntlookup.c fasta_map.a && diff -s ground_truth_seriasm/seeds.rc.after.paf <(./test_xseed.out | sort -h -k1,1 -k6,6)

