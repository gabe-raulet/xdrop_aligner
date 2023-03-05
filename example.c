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

    int err;

    xseed_t xseed;
    xdrop_aligner_t xalign;
    xdrop_score_scheme_t scheme;

    err = xdrop_aligner_set(&xalign, seqQ, seqT);

    if (err == -1) err("xdrop_aligner_set");

    err = xseed_set(&xseed, 5, 3, 3);

    if (err == -1) err("xseed_set");

    err = xseed_check_valid(xseed, xalign);

    if (err == -1) err("xseed_check_valid");

    err = xdrop_score_scheme_set(&scheme, 1, -1, -1, 15);

    if (err == -1) err("xdrop_score_scheme_set");

    err = xdrop_aligner_clear(&xalign);

    if (err == -1) err("xdrop_aligner_clear");

    return 0;
}

