#pragma once

#include <stdint.h>

typedef struct
{
    int8_t *seqQ, *seqT;
    int lenQ, begQ, endQ;
    int lenT, begT, endT;
    int rc;
} xseed_t;
