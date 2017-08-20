
#ifndef _SNOT_PDL_VAX2IEEE_H
#define _SNOT_PDL_VAX2IEEE_H

#include "convert_vax_data.h"


short  vax2s (short  le_short);
int vax2l ( int le_long);
float vax2f (float le_flt);
double vax2d (double le_dbl);


/*
unsigned long int Endian_DWord_Conversion(unsigned long int dword);
unsigned short int Endian_Word_Conversion(unsigned short int word); 
*/

int is_big_endian(void);

void SwapBytes(void *source, int size);


#endif /* _SNOT_PDL_VAX2IEEE_H */
