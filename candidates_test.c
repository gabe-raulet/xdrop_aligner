#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "fasta_map.h"
#include "xdrop_aligner.h"
#include "usage.h"

static inline int maximum(int a, int b) { return a > b? a : b; }

int numreads;
char *buf, **seqs;

int main(int argc, char *argv[])
{
    fasta_map_t fmap = fasta_map_create("reads.fa");
    numreads = loadseqs(fmap, &buf, &seqs);
    fasta_map_free(fmap);

    xdrop_score_scheme_t scheme;
    assert(xdrop_score_scheme_set(&scheme, 1, -1, -1, 15) != -1);

    FILE *f = fopen("candidates.mtx", "r");

    int idxQ, idxT, begQs[2], begTs[2];

    while (fscanf(f, "%d %d %d %d %d %d\n", &idxQ, &idxT, begQs, begTs, begQs+1, begTs+1) > 0)
    {
        if (idxQ == idxT) continue;

        idxQ--, idxT--;

        char *seqQ = seqs[idxQ];
        char *seqT = seqs[idxT];

        int rc;
        xseed_t xseed, result;
        xdrop_seq_pair_t xalign;

        assert(xdrop_seq_pair_set(&xalign, seqQ, seqT) != -1);

        for (int cnt = 0; cnt < 2; ++cnt)
        {
            if (begQs[cnt] == 0 || begTs[cnt] == 0)
                continue;

            if (xseed_set(&xseed, xalign, begQs[cnt], begTs[cnt], 31) == -1)
                continue;

            int score = xdrop_seed_and_extend(xalign, xseed, scheme, &result);

            assert(score != -1);

            int maplen = maximum(result.endT - result.begT, result.endQ - result.begQ);

            printf("%d\t%d\t%d\t%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\n", idxQ, xalign.lenQ, result.begQ, result.endQ, "+-"[xseed.rc], idxT, xalign.lenT, result.begT, result.endT, score, maplen);
        }

        /* int err = xseed_set(&xseed, xalign, begQ2, begT2, 31);                                                                                                                           */

        /* if (err == -1)                                                                                                                                                                   */
        /* {                                                                                                                                                                                */
        /*     printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", idxQ, xalign.lenQ, result.begQ, result.endQ,  idxT, xalign.lenT, result.begT, result.endT);                                       */
        /*     exit(-1);                                                                                                                                                                    */
        /* }                                                                                                                                                                                */

        /* score = xdrop_seed_and_extend(xalign, xseed, scheme, &result);                                                                                                                   */
        /* assert(score != -1);                                                                                                                                                             */

        /* maplen = maximum(result.endT - result.begT, result.endQ - result.begQ);                                                                                                          */

        /* printf("%d\t%d\t%d\t%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\n", idxQ, xalign.lenQ, result.begQ, result.endQ, "+-"[xseed.rc], idxT, xalign.lenT, result.begT, result.endT, score, maplen); */

        assert(xdrop_seq_pair_clear(&xalign) != -1);
    }

    fclose(f);
    free(seqs);
    free(buf);

    return 0;
}
