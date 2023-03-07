#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "fasta_map.h"
#include "xdrop_aligner.h"
#include "ntlookup.h"
#include "usage.h"

#define max(a, b) (((a) > (b)) ? (a) : (b))
#define min(a, b) (((a) < (b)) ? (a) : (b))

/* static inline int max(int a, int b) { return a > b? a : b; } */
/* static inline int min(int a, int b) { return a < b? a : b; } */

int numreads;
char *buf, **seqs;
int *seq_lengths;
int8_t *bufi, **seqsi;

int main(int argc, char *argv[])
{
    fasta_map_t fmap = fasta_map_create("reads.fa");
    numreads = loadseqs(fmap, &buf, &seqs);

    seq_lengths = malloc(numreads * sizeof(int));

    int totbases = 0;

    for (size_t i = 0; i < (size_t)numreads; ++i)
    {
        size_t len;
        fasta_map_query_seq_length_by_id(fmap, i, &len);
        seq_lengths[i] = (int)len;
        totbases += len;
    }

    fasta_map_free(fmap);

    bufi = malloc(totbases);
    seqsi = malloc(numreads * sizeof(int8_t*));

    int8_t *offset = bufi;

    for (size_t i = 0; i < (size_t)numreads; ++i)
    {
        for (int j = 0; j < seq_lengths[i]; ++j)
        {
            offset[j] = NT_LOOKUP_CODE(seqs[i][j]);
        }
        seqsi[i] = offset;
        offset += seq_lengths[i];
    }

    free(seqs);
    free(buf);

    xdrop_score_scheme_t scheme;
    assert(xdrop_score_scheme_set(&scheme, 1, -1, -1, 15) != -1);

    FILE *f = fopen("candidates.mtx", "r");

    int idxQ, idxT, begQs[2], begTs[2];

    while (fscanf(f, "%d %d %d %d %d %d\n", &idxQ, &idxT, begQs, begTs, begQs+1, begTs+1) > 0)
    {
        if (idxQ == idxT) continue;

        idxQ--, idxT--;

        int rc;
        xseed_t xseed, result;
        xdrop_seq_pair_t xalign;

        assert(xdrop_seq_pair_set_ptrs(&xalign, seqsi[idxQ], seqsi[idxT], seq_lengths[idxQ], seq_lengths[idxT]) != -1);

        char classification[128];

        for (int cnt = 0; cnt < 2; ++cnt)
        {
            if (begQs[cnt] == 0 || begTs[cnt] == 0)
                continue;

            if (xseed_set(&xseed, xalign, begQs[cnt], begTs[cnt], 31) == -1)
                continue;

            int score = xdrop_seed_and_extend(xalign, xseed, scheme, &result);

            assert(score != -1);


            int begTr = xseed.rc? xalign.lenT - result.endT : result.begT;
            int endTr = xseed.rc? xalign.lenT - result.begT : result.endT;

            int maplen = max(result.endT - result.begT, result.endQ - result.begQ);
            int overhang  = min(result.begQ, begTr) + min(xalign.lenQ - result.endQ, xalign.lenT - result.endT);

            if (overhang > min(1000, maplen * 0.8))
            {
                sprintf(classification, "internal_match");
            }
            else if (result.begQ <= begTr && xalign.lenQ - result.endQ <= xalign.lenT - endTr)
            {
                sprintf(classification, "first_contained");
            }
            else if (result.begQ >= begTr && xalign.lenQ - result.endQ >= xalign.lenT - endTr)
            {
                sprintf(classification, "second_contained");
            }
            else if (result.begQ > result.begT)
            {
                sprintf(classification, "first_to_second_overlap");
            }
            else
            {
                sprintf(classification, "second_to_first_overlap");
            }

            printf("%d\t%d\t%d\t%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n", idxQ, xalign.lenQ, result.begQ, result.endQ, "+-"[xseed.rc], idxT, xalign.lenT, result.begT, result.endT, score, maplen, classification);
        }

        assert(xdrop_seq_pair_clear_ptrs(&xalign) != -1);
    }

    fclose(f);
    free(seqsi);
    free(bufi);

    return 0;
}
