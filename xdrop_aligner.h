#pragma once

#include <stdint.h>

/*
 * xdrop sequence pair structure and functions
 */

typedef struct
{
    int8_t *seqQ, *seqT, *seqTr;
    int lenQ, lenT;
} xdrop_seq_pair_t;

int xdrop_seq_pair_set(xdrop_seq_pair_t *xalign, char const *seqQ, char const *seqT);
int xdrop_seq_pair_clear(xdrop_seq_pair_t *xalign);

/*
 * seed structure and functions
 */

typedef struct
{
    int begQ, endQ;
    int begT, endT;
    int rc;
} xseed_t;

int xseed_set(xseed_t *xseed, xdrop_seq_pair_t const refpair, int begQ, int begT, int seedlen);
int xseed_seedlen(xseed_t const xseed);

/*
 * scoring structure and functions
 */

typedef struct
{
    int mat;
    int mis;
    int gap;
    int dropoff;
} xdrop_score_scheme_t;

int xdrop_score_scheme_set(xdrop_score_scheme_t *scheme, int mat, int mis, int gap, int dropoff);

/*
 * alignment functions
 */

int
xdrop_seed_and_extend_l
(
    xdrop_seq_pair_t const xalign,
    xdrop_score_scheme_t const scheme,
    xseed_t xseed,
    int *begQ_ext,
    int *begT_ext
);

int
xdrop_seed_and_extend_r
(
    xdrop_seq_pair_t const xalign,
    xdrop_score_scheme_t const scheme,
    xseed_t xseed,
    int *endQ_ext,
    int *endT_ext
);

int
xdrop_seed_and_extend
(
    xdrop_seq_pair_t const xalign,
    xseed_t const xseed,
    xdrop_score_scheme_t const scheme,
    xseed_t *result
);
