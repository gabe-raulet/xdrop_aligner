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
/* } xdrop_seq_pair_t;       */

int xdrop_seq_pair_set(xdrop_seq_pair_t *xalign, char const *seqQ, char const *seqT)
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

int xdrop_seq_pair_clear(xdrop_seq_pair_t *xalign)
{
    if (!xalign) return -1;

    if (xalign->seqQ) free(xalign->seqQ);
    if (xalign->seqT) free(xalign->seqT);
    if (xalign->seqTr) free(xalign->seqTr);

    memset(xalign, 0, sizeof(xdrop_seq_pair_t));

    return 0;
}

/* typedef struct      */
/* {                   */
/*     int begQ, endQ; */
/*     int begT, endT; */
/* } xseed_t;          */

int xseed_set(xseed_t *xseed, xdrop_seq_pair_t const refpair, int begQ, int begT, int seedlen)
{
    /*
     * Seedlen must be odd, otherwise we might have reverse complement palindromes.
     */
    if (!xseed || !(seedlen&1)) return -1;

    /*
     * Check that the query sequence coordinates are within the logical bounds.
     */
    if (begQ < 0 || begQ + seedlen >= refpair.lenQ)
        return -1;

    /*
     * Check that the target sequence coordinates are within the logical bounds.
     */
    if (begT < 0 || begT + seedlen >= refpair.lenT)
        return -1;

    /*
     * Assuming that the seed coordinates are correct (which we check just after this),
     * the pairs can be determined to be on opposite strands by simply looking at the
     * nucleotide right in the middle of the seed on each sequence and seeing whether
     * the are Watson-Crick complements or not.
     */
    int rc = (refpair.seqQ[begQ + (seedlen>>1)] == 3 - refpair.seqT[begT + (seedlen>>1)]);

    for (int i = 0; i < seedlen; ++i)
    {
        if (refpair.seqQ[begQ + i] != (rc? 3 - refpair.seqT[begT + seedlen - 1 - i] : refpair.seqT[begT + i]))
            return -1;
    }

    xseed->begQ = begQ;
    xseed->endQ = xseed->begQ + seedlen;

    xseed->begT = rc? refpair.lenT - begT - seedlen : begT;
    xseed->endT = xseed->begT + seedlen;

    xseed->rc = rc;

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

int extend_seed_one_direction(const int8_t *target_seq, const int8_t *query_seq, int target_len, int query_len, xseed_t *xseed, int extleft, int mat, int mis, int gap, int xdrop)
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
    int target_offset = xseed->endT;
    int query_offset = xseed->endQ;

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

    if (best_ext_score != undef)
    {
        if (extleft)
        {
            xseed->begT -= best_ext_row;
            xseed->begQ -= best_ext_col;
        }
        else
        {
            xseed->endT += best_ext_row;
            xseed->endQ += best_ext_col;
        }
    }

    vec_free(ad1);
    vec_free(ad2);
    vec_free(ad3);

    return best_ext_score;
}


int xdrop_seed_and_extend(xdrop_seq_pair_t const xalign, xseed_t const xseed, xdrop_score_scheme_t const scheme, xseed_t *result)
{
    int seedlen = xseed_seedlen(xseed);

    if (seedlen == -1)
        return -1;

    int begQ_ext, begT_ext, lscore;
    int endQ_ext, endT_ext, rscore;

    lscore = xdrop_seed_and_extend_l(xalign, scheme, xseed, &begQ_ext, &begT_ext);
    rscore = xdrop_seed_and_extend_r(xalign, scheme, xseed, &endQ_ext, &endT_ext);

    int score = lscore + rscore + scheme.mat * seedlen;

    result->begQ = begQ_ext;
    result->endQ = endQ_ext;

    result->begT = xseed.rc? xalign.lenT - endT_ext : begT_ext;
    result->endT = xseed.rc? xalign.lenT - begT_ext : endT_ext;

    return score;
}

int xdrop_seed_and_extend_l(xdrop_seq_pair_t const xalign, xdrop_score_scheme_t const scheme, xseed_t xseed, int *begQ_ext, int *begT_ext)
{
    int8_t const *seqT = xseed.rc? xalign.seqTr : xalign.seqT;

    int lscore = extend_seed_one_direction(seqT, xalign.seqQ, xseed.begT, xseed.begQ, &xseed, 1, scheme.mat, scheme.mis, scheme.gap, scheme.dropoff);

    *begQ_ext = xseed.begQ;
    *begT_ext = xseed.begT;

    return lscore;
}

int xdrop_seed_and_extend_r(xdrop_seq_pair_t const xalign, xdrop_score_scheme_t const scheme, xseed_t xseed, int *endQ_ext, int *endT_ext)
{
    int8_t const *seqT = xseed.rc? xalign.seqTr : xalign.seqT;

    int rscore = extend_seed_one_direction(seqT, xalign.seqQ, xalign.lenT - xseed.endT, xalign.lenQ - xseed.endQ, &xseed, 0, scheme.mat, scheme.mis, scheme.gap, scheme.dropoff);

    *endQ_ext = xseed.endQ;
    *endT_ext = xseed.endT;

    return rscore;
}
