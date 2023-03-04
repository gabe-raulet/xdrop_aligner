
# for reverse complement

base_for = "ACGT"
base_rev = "TGCA"
comp_tab = str.maketrans(base_for, base_rev)

def revcomp(s): return s.translate(comp_tab)[::-1]

# for reading mixed upper/lower/'N' nucleotide sequences

char_in  = "acgtnACGTN"
char_out = "ACGTAACGTA"
nt_tab = str.maketrans(char_in, char_out)

def nt_normalize(s): return s.translate(nt_tab)

# examples

seqQ = "agtggCAAga"
seqT =   "tggCAAgaccata"
twin = "tatggtcTTGcca"
