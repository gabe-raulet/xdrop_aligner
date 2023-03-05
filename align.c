#include "align.h"
#include "vec.h"
#include "usage.h"
#include <limits.h>
#include <stdlib.h>
#include <string.h>

#define max(a, b) (((a) < (b))? (b) : (a))
#define min(a, b) (((a) > (b))? (b) : (a))

int extend_seed_one_direction(const char *target_seq, const char *query_seq, int target_len, int query_len, ext_seed_t *xseed, int extleft, int mat, int mis, int gap, int xdrop)
{
    /* printf("target_len=%d, query_len=%d, extleft=%d, tbeg=%d, tend=%d, qbeg=%d, qend=%d\n", */
            /* target_len, query_len, extleft, xseed->tbeg, xseed->tend, xseed->qbeg, xseed->qend); */

    int cols = query_len + 1;
    int rows = target_len + 1;

    /* printf("cols=%d, rows=%d\n", cols, rows); */

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

    /* printf("target_offset=%d, query_offset=%d\n", target_offset, query_offset); */

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

    /* printf("best_ext_score=%d, tbeg_best=%d, tend_best=%d, qbeg_best=%d, qend_best=%d\n", best_ext_score, xseed->tbeg_best, xseed->tend_best, xseed->qbeg_best, xseed->qend_best); */

    return best_ext_score;
}

int extend_seed(const char *target_seq, const char *query_seq, int target_len, int query_len, ext_seed_t *xseed, int mat, int mis, int gap, int xdrop)
{
    /* char *target_prefix = strndup(target_seq, xseed->tbeg); */
    /* char *query_prefix = strndup(query_seq, xseed->qbeg); */

    int lscore = extend_seed_one_direction(target_seq, query_seq, xseed->tbeg, xseed->qbeg, xseed, 1, mat, mis, gap, xdrop);

    /* free(target_prefix); */
    /* free(query_prefix); */

    /* char *target_suffix = strndup(target_seq + xseed->tend, target_len - xseed->tend); */
    /* char *query_suffix = strndup(query_seq + xseed->qend, query_len - xseed->qend); */

    int rscore = extend_seed_one_direction(target_seq, query_seq, target_len - xseed->tend, query_len - xseed->qend, xseed, 0, mat, mis, gap, xdrop);

    /* free(target_suffix); */
    /* free(query_suffix); */

    return lscore + rscore;
}
