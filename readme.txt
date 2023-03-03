/*------- XSeed --------*/

The purpose of this module is to develop a very painless method to run X-drop seed-
and-extend alignments for nucleotide sequences, in the C programming language. A
seed-and-extend alignment takes the following as sufficient inputs:

    1. Two nucleotide sequences. We name them seqQ and seqT, where the Q/T references
       the common names provided when two sequences are aligned: seqQ is the "query"
       sequence and seqT is the "target" sequence. It is assumed that when aligning
       two sequences that originate from opposite strands, it is always the target
       sequence which is reverse complemented. Hence, a same-strand alignment would
       be between an infix of seqQ and an infix of seqT, whereas an opposite-strand
       alignment would be between an infix of seqQ and an infix of rc(seqT), where
       for any nucleotide sequences A, rc(A) denotes the reverse complement of A.

    2. In addition to the two sequences, we need the coordinates on each sequence of
       the shared (common) seed between them. This is fundamental to any seed-and-extend
       alignment algorithm. Seed-and-extend methods almost always rely upon a pipeline
       which searches for common seeds (infix subsequences) that are shared between
       different sequences. The idea is that if there is a common infix seed shared between
       two distinct sequences, then there is a possibility that the neighboring sequences
       on the left and right flanks of each seed (in their respective sequences) have
       a high degree of sequence similarity.

       The coordinates can be supplied two ways, although the first is the most logical:

          i) begQ, begT, and seedlen. Given these coordinates, it is implied that
             seqQ[begQ..begQ+seedlen-1] $= seqT[begT..begT+seedlen-1], where I am
             using the s

