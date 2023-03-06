#!/bin/bash

clang -g -O0 -o test_xseed.out test_xseed.c xdrop_aligner.c ntlookup.c fasta_map.a

