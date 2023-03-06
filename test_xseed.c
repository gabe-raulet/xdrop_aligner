#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "fasta_map.h"
#include "xdrop_aligner.h"
#include "usage.h"

char const *fasta_fname = "ground_truth_seriasm/reads.fa";
char const *seeds_fname = "ground_truth_seriasm/seeds.sorted.before.paf";

int numreads;
char *buf, **seqs;

int main(int argc, char *argv[])
{
    fasta_map_t fmap = fasta_map_create(fasta_fname);
    numreads = loadseqs(fmap, &buf, &seqs);
    fasta_map_free(fmap);

    FILE *f = fopen(seeds_fname, "r");

    int idxQ, lenQ, begQ, endQ, idxT, lenT, begT, endT;

    while (fscanf(f, "%d %d %d %d %*c %d %d %d %d %*d %*d %*d %*s\n", &idxQ, &lenQ, &begQ, &endQ, &idxT, &lenT, &begT, &endT) > 0)
    {
        char *seqQ = seqs[idxQ-1];
        char *seqT = seqs[idxT-1];
        int lenQ_chk = strlen(seqQ);
        int lenT_chk = strlen(seqT);

        assert(lenQ_chk == lenQ);
        assert(lenT_chk == lenT);
    }

    fclose(f);
    free(buf);
    free(seqs);

    return 0;
}
