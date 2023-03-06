#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "xdrop_aligner.h"
#include "usage.h"

int main(int argc, char *argv[])
{
    char const *seqQ = "agtggCAAga";
    char const *seqT =   "tggCAAgaccata";
    char const *twin = "tatgggtgTTGcca";

    xseed_t xseed;
    xdrop_aligner_t xalign;
    xdrop_score_scheme_t scheme;

    assert(xdrop_aligner_set(&xalign, seqQ, seqT) != -1);
    assert(xseed_set(&xseed, xalign, 5, 3, 3, 0) != -1);
    assert(xseed_check_valid(xseed, xalign) != -1);
    assert(xdrop_score_scheme_set(&scheme, 1, -1, -1, 15) != -1);

    xseed_t result;

    printf("%d\t%d\t%d\t%d\t%d\t%d\n", xalign.lenQ, xseed.begQ, xseed.endQ, xalign.lenT, xseed.begT, xseed.endT);

    int score = xdrop_seed_and_extend(xalign, xseed, scheme, &result);

    printf("%d\t%d\t%d\t%d\t%d\t%d\n", xalign.lenQ, result.begQ, result.endQ, xalign.lenT, result.begT, result.endT);

    assert(xdrop_aligner_clear(&xalign) != -1);

    return 0;
}

