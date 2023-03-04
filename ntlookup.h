#pragma once

extern char const nt_lookup_char[4];
extern char const nt_lookup_code[256];

#define NT_LOOKUP_CHAR(c) (nt_lookup_char[(int)(c)])
#define NT_LOOKUP_CODE(c) (nt_lookup_code[(int)(c)])
