#include "align.h"
#include "vec.h"
#include "seq.h"
#include "usage.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

static inline int max(int a, int b) { return a > b? a : b; }
static inline int min(int a, int b) { return a < b? a : b; }

int main(int argc, char *argv[])
{
    char const *seqQ = "agtggCAAga";
    char const *seqT =   "tggCAAgaccata";
    char const *twin = "tatgggtgTTGcca";

    int lenQ = strlen(seqQ);
    int lenT = strlen(seqT);

    ext_seed_t xseed;
    xseed.qbeg = xseed.qbeg_best = 5;
    xseed.tbeg = xseed.tbeg_best = 3;
    xseed.qend = xseed.qend_best = xseed.qbeg + 3;
    xseed.tend = xseed.tend_best = xseed.tbeg + 3;

    printf("%d\t%d\t%d\t%d\t%d\t%d\n", lenQ, xseed.qbeg_best, xseed.qend_best, lenT, xseed.tbeg_best, xseed.tend_best);

    int score = extend_seed(seqT, seqQ, lenT, lenQ, &xseed, 1, -1, -1, 15);

    printf("%d\t%d\t%d\t%d\t%d\t%d\n", lenQ, xseed.qbeg_best, xseed.qend_best, lenT, xseed.tbeg_best, xseed.tend_best);

    return 0;
}
