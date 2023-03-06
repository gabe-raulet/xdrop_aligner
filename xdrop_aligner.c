#include "xdrop_aligner.h"
#include "ntlookup.h"
#include "vec.h"
#include "usage.h"
#include <stdlib.h>
#include <string.h>
#include <limits.h>

static inline int max(int a, int b) { return a > b? a : b; }
static inline int min(int a, int b) { return a < b? a : b; }

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
    xalign->seqTr = malloc(xalign->lenT);

    int i;

    for (i = 0; i < xalign->lenQ; ++i)
    {
        xalign->seqQ[i] = NT_LOOKUP_CODE(seqQ[i]);
    }

    for (i = 0; i < xalign->lenT; ++i)
    {
        xalign->seqT[i] = NT_LOOKUP_CODE(seqT[i]);
        xalign->seqTr[xalign->lenT - 1 - i] = 3 - xalign->seqT[i];
    }

    return 0;
}

int xdrop_aligner_clear(xdrop_aligner_t *xalign)
{
    if (!xalign) return -1;

    if (xalign->seqQ) free(xalign->seqQ);
    if (xalign->seqT) free(xalign->seqT);
    if (xalign->seqTr) free(xalign->seqTr);

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

typedef struct
{
    int tbeg, qbeg, tend, qend;
    int tbeg_best, qbeg_best, tend_best, qend_best;
} ext_seed_t;

int extend_seed_one_direction(const int8_t *target_seq, const int8_t *query_seq, int target_len, int query_len, ext_seed_t *xseed, int extleft, int mat, int mis, int gap, int xdrop)
{
    int cols = query_len + 1;
    int rows = target_len + 1;

    if (rows == 1 || cols == 1) return 0;

    int len = 2 * max(cols, rows);
    int min_err_score = INT_MIN / len;
    gap = max(gap, min_err_score);
    mis = max(mis, min_err_score);
    int undef = INT_MIN - gap;

    vec_t(int) ad1, ad2, ad3, tmp;
    vec_init(ad1);
    vec_init(ad2);
    vec_init(ad3);

    int min_col = 1, max_col = 2;
    int offset1 = 0, offset2 = 0, offset3 = 0;

    vec_reserve(ad2, 2);
    vec_at(ad2, 0) = 0;
    ad2.l = 1;

    int best_ext_col = 0, best_ext_row = 0, best_ext_score = 0;

    vec_reserve(ad3, 4);
    vec_at(ad3, 0) = vec_at(ad3, 1) = (-gap > xdrop)? undef : gap;
    ad3.l = 2;

    int ad_no = 1, best = 0;
    int target_offset = xseed->tend;
    int query_offset = xseed->qend;

    while (min_col < max_col)
    {
        ++ad_no;
        tmp = ad1;
        ad1 = ad2;
        ad2 = ad3;
        ad3 = tmp;

        offset1 = offset2;
        offset2 = offset3;
        offset3 = min_col - 1;

        vec_reserve(ad3, max_col+1-offset3);
        vec_at(ad3, 0) = vec_at(ad3, max_col - offset3) = undef;
        ad3.l = max_col+1-offset3;

        if (ad_no * gap > best - xdrop)
        {
            if (offset3 == 0) vec_at(ad3, 0) = ad_no * gap;
            if (ad_no - max_col == 0) vec_at(ad3, max_col - offset3) = ad_no * gap;
        }

        int ad_best = ad_no * gap;

        for (int col = min_col; col < max_col; ++col)
        {
            int i3 = col - offset3;
            int i2 = col - offset2;
            int i1 = col - offset1;

            int query_pos, target_pos;

            query_pos = extleft? cols - 1 - col : col - 1 + query_offset;
            target_pos = extleft? rows - 1 + col - ad_no : ad_no - col - 1 + target_offset;

            int temp = max(ad2.a[i2-1], ad2.a[i2]) + gap;
            int temp2 = vec_at(ad1, i1-1) + (query_seq[query_pos] == target_seq[target_pos]? mat : mis);
            temp = max(temp, temp2);

            if (temp < best - xdrop)
            {
                vec_at(ad3, i3) = undef;
            }
            else
            {
                vec_at(ad3, i3) = temp;
                ad_best = max(ad_best, temp);
            }

            if (temp > best)
            {
                best_ext_col = col;
                best_ext_row = ad_no - best_ext_col;
                best_ext_score = vec_at(ad3, best_ext_col - offset3);

                if (best_ext_score != temp) bug("assert");
            }
        }

        best = max(best, ad_best);

        while (min_col - offset3 < vec_size(ad3) && vec_at(ad3, min_col - offset3) == undef &&
               min_col - offset2 - 1 < vec_size(ad2) && vec_at(ad2, min_col - offset2 - 1) == undef)
        {
            ++min_col;
        }

        while (max_col - offset3 > 0 && vec_at(ad3, max_col - offset3 - 1) == undef && vec_at(ad2, max_col - offset2 - 1) == undef)
            --max_col;

        ++max_col;

        min_col = max(min_col, ad_no + 2 - rows);
        max_col = min(max_col, cols);
    }

    int ext_col = vec_size(ad3) + offset3 - 2;
    int ext_row = ad_no - ext_col;
    int ext_score = vec_at(ad3, ext_col - offset3);

    if (ext_score == undef)
    {
        if (vec_at(ad2, ad2.l-2) != undef)
        {
            ext_col = vec_size(ad2) + offset2 - 2;
            ext_row = ad_no - 1 - ext_col;
            ext_score = vec_at(ad2, ext_col - offset2);
        }
        else if (vec_size(ad2) > 2 && vec_at(ad2, ad2.l-3) != undef)
        {
            ext_col = vec_size(ad2) + offset2 - 3;
            ext_row = ad_no - 1 - ext_col;
            ext_score = vec_at(ad2, ext_col - offset3);
        }
    }

    if (ext_score == undef)
    {
        for (int i = 0; i < vec_size(ad1); ++i)
        {
            if (vec_at(ad1, i) > ext_score)
            {
                ext_score = vec_at(ad1, i);
                ext_col = i + offset1;
                ext_row = ad_no - 2 - ext_col;
            }
        }
    }

    if (ext_score != undef)
    {
        if (extleft)
        {
            xseed->tbeg -= ext_row;
            xseed->qbeg -= ext_col;
        }
        else
        {
            xseed->tend += ext_row;
            xseed->qend += ext_col;
        }
    }

    if (best_ext_score != undef)
    {
        if (extleft)
        {
            xseed->tbeg_best -= best_ext_row;
            xseed->qbeg_best -= best_ext_col;
        }
        else
        {
            xseed->tend_best += best_ext_row;
            xseed->qend_best += best_ext_col;
        }
    }

    vec_free(ad1);
    vec_free(ad2);
    vec_free(ad3);

    return best_ext_score;
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

    result->begQ = begQ_ext;
    result->endQ = endQ_ext;
    result->begT = begT_ext;
    result->endT = endT_ext;

    return score;
}

int xdrop_seed_and_extend_l(xdrop_aligner_t const xalign, xdrop_score_scheme_t const scheme, int begQ, int begT, int endQ, int endT, int *begQ_ext, int *begT_ext)
{
    ext_seed_t xseed = {begT, begQ, endT, endQ, begT, begQ, endT, endQ};

    int lscore = extend_seed_one_direction(xalign.seqT, xalign.seqQ, begT, begQ, &xseed, 1, scheme.mat, scheme.mis, scheme.gap, scheme.dropoff);

    *begQ_ext = xseed.qbeg_best;
    *begT_ext = xseed.tbeg_best;

    return lscore;
}

int xdrop_seed_and_extend_r(xdrop_aligner_t const xalign, xdrop_score_scheme_t const scheme, int begQ, int begT, int endQ, int endT, int *endQ_ext, int *endT_ext)
{
    ext_seed_t xseed = {begT, begQ, endT, endQ, begT, begQ, endT, endQ};

    int rscore = extend_seed_one_direction(xalign.seqT, xalign.seqQ, xalign.lenT - endT, xalign.lenQ - endQ, &xseed, 0, scheme.mat, scheme.mis, scheme.gap, scheme.dropoff);

    *endQ_ext = xseed.qend_best;
    *endT_ext = xseed.tend_best;

    return rscore;
}
