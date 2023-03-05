#include "xdrop_aligner.h"
#include "ntlookup.h"
#include <stdlib.h>
#include <string.h>

/* typedef struct           */
/* {                        */
/*     int8_t *seqQ, *seqT; */
/*     int lenQ, lenT;      */
/* } xdrop_aligner_t;       */

int xdrop_aligner_set(xdrop_aligner_t *xalign, char const *seqQ, char const *seqT)
{
    if (!xalign || !seqQ || !seqT)
        return -1;

    xalign->lenQ = strlen(seqQ);
    xalign->lenT = strlen(seqT);

    xalign->seqQ = malloc(xalign->lenQ);
    xalign->seqT = malloc(xalign->lenT);

    int i;

    for (i = 0; i < xalign->lenQ; ++i)
    {
        xalign->seqQ[i] = NT_LOOKUP_CODE(seqQ[i]);
    }

    for (i = 0; i < xalign->lenT; ++i)
    {
        xalign->seqT[i] = NT_LOOKUP_CODE(seqT[i]);
    }

    return 0;
}

int xdrop_aligner_clear(xdrop_aligner_t *xalign)
{
    if (!xalign) return -1;

    if (xalign->seqQ) free(xalign->seqQ);
    if (xalign->seqT) free(xalign->seqT);

    memset(xalign, 0, sizeof(xdrop_aligner_t));

    return 0;
}

/* typedef struct      */
/* {                   */
/*     int begQ, endQ; */
/*     int begT, endT; */
/* } xseed_t;          */

int xseed_set(xseed_t *xseed, int begQ, int begT, int seedlen)
{
    if (!xseed) return -1;

    xseed->begQ = begQ;
    xseed->begT = begT;
    xseed->endQ = begQ + seedlen;
    xseed->endT = begT + seedlen;

    return 0;
}

int xseed_check_valid(xseed_t const xseed, xdrop_aligner_t const xalign)
{
    /*
     * Check that the xseed coordinates are consistent with the sequences
     * in the xalign structure.
     */

    if (xseed.begQ < 0 || xseed.endQ >= xalign.lenQ)
        return -1;

    if (xseed.begT < 0 || xseed.endT >= xalign.lenT)
        return -1;

    /*
     * Check that the xseed coordinates correspond to equal length seeds
     * for both the query and the target.
     */

    int seedlen = xseed_seedlen(xseed);

    if (seedlen == -1)
        return -1;

    /*
     * Check that the seed sequences are identical on the query and target
     * sequences.
     */

    for (int i = 0; i < seedlen; ++i)
    {
        if (xalign.seqQ[xseed.begQ + i] != xalign.seqT[xseed.begT + i])
            return -1;
    }

    return 0;
}

int xseed_seedlen(xseed_t const xseed)
{
    int seedlenQ = xseed.endQ - xseed.begQ;
    int seedlenT = xseed.endT - xseed.begT;

    return seedlenQ != seedlenT? -1 : seedlenQ;
}

/* typedef struct          */
/* {                       */
/*     int mat;            */
/*     int mis;            */
/*     int gap;            */
/*     int dropoff;        */
/* } xdrop_score_scheme_t; */

int xdrop_score_scheme_set(xdrop_score_scheme_t *scheme, int mat, int mis, int gap, int dropoff)
{
    if (!scheme) return -1;

    scheme->mat = mat;
    scheme->mis = mis;
    scheme->gap = gap;
    scheme->dropoff = dropoff;

    return 0;
}

int xdrop_seed_and_extend(xdrop_aligner_t const xalign, xseed_t const xseed, xdrop_score_scheme_t const scheme, xseed_t *result)
{
    int seedlen = xseed_seedlen(xseed);

    if (seedlen == -1)
        return -1;

    int begQ, begT, begQ_ext, begT_ext, lscore;
    int endQ, endT, endQ_ext, endT_ext, rscore;

    begQ = xseed.begQ;
    begT = xseed.begT;
    endQ = xseed.endQ;
    endT = xseed.endT;

    lscore = xdrop_seed_and_extend_l(xalign, scheme, begQ, begT, endQ, endT, &begQ_ext, &begT_ext);
    rscore = xdrop_seed_and_extend_r(xalign, scheme, begQ, begT, endQ, endT, &endQ_ext, &endT_ext);

    int score = lscore + rscore + scheme.mat * seedlen;

    return score;
}

int xdrop_seed_and_extend_l(xdrop_aligner_t const xalign, xdrop_score_scheme_t const scheme, int begQ, int begT, int endQ, int endT, int *begQ_ext, int *begT_ext);
int xdrop_seed_and_extend_r(xdrop_aligner_t const xalign, xdrop_score_scheme_t const scheme, int begQ, int begT, int endQ, int endT, int *endQ_ext, int *endT_ext);
