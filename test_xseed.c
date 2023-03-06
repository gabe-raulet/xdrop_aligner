#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "fasta_map.h"
#include "xdrop_aligner.h"
#include "usage.h"

static inline int maximum(int a, int b) { return a > b? a : b; }

char const *fasta_fname = "ground_truth_seriasm/reads.fa";
char const *seeds_fname = "ground_truth_seriasm/seeds.before.paf";

int numreads;
char *buf, **seqs;

int main(int argc, char *argv[])
{
    fasta_map_t fmap = fasta_map_create(fasta_fname);
    numreads = loadseqs(fmap, &buf, &seqs);
    fasta_map_free(fmap);

    xdrop_score_scheme_t scheme;
    assert(xdrop_score_scheme_set(&scheme, 1, -1, -1, 15) != -1);

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

        xseed_t xseed;
        xdrop_seq_pair_t xalign;
        int rc;

        assert(xdrop_seq_pair_set(&xalign, seqQ, seqT) != -1);

        assert(xseed_set(&xseed, xalign, begQ, begT, 31) != -1);

        xseed_t result;

        int score = xdrop_seed_and_extend(xalign, xseed, scheme, &result);

        assert(score != -1);

        int maplen = maximum(result.endT - result.begT, result.endQ - result.begQ);

        printf("%d\t%d\t%d\t%d\t+\t%d\t%d\t%d\t%d\t%d\t%d\t255\tAFTER\n", idxQ, xalign.lenQ, result.begQ, result.endQ,
                                                                          idxT, xalign.lenT, result.begT, result.endT, score, maplen);

        assert(xdrop_seq_pair_clear(&xalign) != -1);

    }

    fclose(f);
    free(buf);
    free(seqs);

    return 0;
}
