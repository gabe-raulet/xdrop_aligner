#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#include <limits.h>
#include "fasta_map.h"
#include "xdrop_aligner.h"
#include "ntlookup.h"
#include "usage.h"

#define max(a, b) (((a) > (b))? (a) : (b))
#define min(a, b) (((a) < (b))? (a) : (b))

int mat = 1;
int mis = -1;
int gap = -1;
int dropoff = 15;

int numreads;
char *buf, **seqs;
int8_t *bufi, **seqsi;
int *seq_lengths;

int usage(char const *prgname)
{
    fprintf(stderr, "Usage: %s [options] <in.fa> <seeds.mtx>\n", prgname);
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -A INT   matching score [%d]\n", mat);
    fprintf(stderr, "    -B INT   mismatch penalty [%d]\n", -mis);
    fprintf(stderr, "    -G INT   gap penalty [%d]\n", -gap);
    fprintf(stderr, "    -x INT   x-drop dropoff [%d]\n", dropoff);
    return -1;
}

int main(int argc, char *argv[])
{
    int c;

    while ((c = getopt(argc, argv, "A:B:G:x:h")) >= 0)
    {
        if      (c == 'A') mat = atoi(optarg);
        else if (c == 'B') mis = -atoi(optarg);
        else if (c == 'G') gap = -atoi(optarg);
        else if (c == 'x') dropoff = atoi(optarg);
    }

    if (optind + 2 > argc) return usage(argv[0]);

    fasta_map_t fmap = fasta_map_create(argv[optind]);
    numreads = loadseqs(fmap, &buf, &seqs);

    size_t totbases;
    double avglen;

    fasta_map_stats(fmap, NULL, &totbases, NULL, NULL, &avglen);

    seq_lengths = malloc(numreads * sizeof(int));

    for (size_t i = 0; i < (size_t)numreads; ++i)
    {
        size_t len;
        fasta_map_query_seq_length_by_id(fmap, i, &len);
        seq_lengths[i] = (int)len;
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

    fprintf(stderr, "Read %d sequences of average length %.3f from FASTA file '%s'\n", numreads, avglen, argv[optind]);

    free(seqs);
    free(buf);

    xdrop_score_scheme_t scheme;
    xdrop_score_scheme_set(&scheme, mat, mis, gap, dropoff);

    FILE *f = fopen(argv[optind+1], "r");

    int idxQ, idxT, begQs[2], begTs[2];

    char line[1024];

    fgets(line, 1024, f);
    fgets(line, 1024, f);

    while (fgets(line, 1024, f) != NULL)
    {
        sscanf(line, "%d %d %d %d %d %d\n", &idxQ, &idxT, begQs, begTs, begQs+1, begTs+1);

        if (idxQ == idxT) continue;

        idxQ--, idxT--;

        int rc;
        xseed_t xseed, result, best_result;
        xdrop_seq_pair_t xalign;

        assert(xdrop_seq_pair_set_ptrs(&xalign, seqsi[idxQ], seqsi[idxT], seq_lengths[idxQ], seq_lengths[idxT]) != -1);

        char classification[128];

        int best_score = INT_MIN;
        int maplen, overhang, begTr, endTr;

        for (int cnt = 0; cnt < 2; ++cnt)
        {
            if (begQs[cnt] == 0 || begTs[cnt] == 0)
                continue;

            if (xseed_set(&xseed, xalign, begQs[cnt], begTs[cnt], 31) == -1)
                continue;

            int score = xdrop_seed_and_extend(xalign, xseed, scheme, &result);

            assert(score != -1);

            if (score < best_score)
                continue;

            best_score = score;

            begTr = xseed.rc? xalign.lenT - result.endT : result.begT;
            endTr = xseed.rc? xalign.lenT - result.begT : result.endT;

            maplen = max(result.endT - result.begT, result.endQ - result.begQ);
            overhang  = min(result.begQ, begTr) + min(xalign.lenQ - result.endQ, xalign.lenT - endTr);

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

            best_result = result;
        }

        printf("%d\t%d\t%d\t%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n", idxQ, xalign.lenQ, best_result.begQ, best_result.endQ, "+-"[xseed.rc], idxT, xalign.lenT, best_result.begT, best_result.endT, best_score, maplen, classification);

        assert(xdrop_seq_pair_clear_ptrs(&xalign) != -1);
    }

    free(seqsi);
    free(bufi);
    free(seq_lengths);

    return 0;
}
