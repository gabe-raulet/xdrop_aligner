#pragma once

typedef struct
{
    int tbeg, qbeg, tend, qend;
    int tbeg_best, qbeg_best, tend_best, qend_best;
} ext_seed_t;

int extend_seed
(
    const char *target_seq,
    const char *query_seq,
    int target_len,
    int query_len,
    ext_seed_t *xseed,
    int mat,
    int mis,
    int gap,
    int xdrop
);

int extend_seed_one_direction
(
    const char *target_seq,
    const char *query_seq,
    int target_len,
    int query_len,
    ext_seed_t *xseed,
    int extleft,
    int mat,
    int mis,
    int gap,
    int xdrop
);
