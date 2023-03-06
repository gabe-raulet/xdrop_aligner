#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "xdrop_aligner.h"
#include "usage.h"

int main(int argc, char *argv[])
{
    char const *seqQ = "agtggCAAga";
    char const *seqT =   "tggCAAgaccata";
    char const *twin = "tatggtcTTGcca";

    xseed_t xseed;
    xdrop_seq_pair_t xalign;
    xdrop_score_scheme_t scheme;

    assert(xdrop_seq_pair_set(&xalign, seqQ, twin) != -1);
    assert(xseed_set(&xseed, xalign, 5, 7, 3) != -1);
    /* assert(xseed_check_valid(xseed, xalign) != -1); */
    assert(xdrop_score_scheme_set(&scheme, 1, -1, -1, 15) != -1);

    xseed_t result;

    printf("%d\t%d\t%d\t%d\t%d\t%d\n", xalign.lenQ, xseed.begQ, xseed.endQ, xalign.lenT, xseed.begT, xseed.endT);

    int score = xdrop_seed_and_extend(xalign, xseed, scheme, &result);

    printf("%d\t%d\t%d\t%d\t%d\t%d\n", xalign.lenQ, result.begQ, result.endQ, xalign.lenT, result.begT, result.endT);

    assert(xdrop_seq_pair_clear(&xalign) != -1);

    return 0;
}

