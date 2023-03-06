#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fasta_map.h"
#include "xdrop_aligner.h"
#include "usage.h"

int numreads;
char *buf, **seqs;

int main(int argc, char *argv[])
{
    fasta_map_t fmap = fasta_map_create("ground_truth_seriasm/reads.fa");
    numreads = loadseqs(fmap, &buf, &seqs);
    fasta_map_free(fmap);

    for (int i = 0; i < numreads; ++i)
    {
        printf("%s\n", seqs[i]);
    }

    free(buf);
    free(seqs);

    return 0;
}
