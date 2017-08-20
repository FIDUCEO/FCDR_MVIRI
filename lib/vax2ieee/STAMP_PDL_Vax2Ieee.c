/*******************************************************************************
*            Copyright Eumetsat 2014
*
*    Filename        :   	SNOT_PDL_Vax2Ieee.c
*
*    Version         :    
*
*    Date            :    
*
*    Status          :    
*
*    Author          :  A. Lattanzio      
*
*    Description     :       
*        
*        
*                                                                            *
*******************************************************************************/

#include "STAMP_PDL_Vax2Ieee.h"


/* The following funxtions are needed to swap from VAX to IEEE in the MFG reader */
short  vax2s (short  le_short)
{
    return(le_short);
}
 int vax2l ( int le_long)
{
    return(le_long);
}

float vax2f( float f )
{
  int n =1;
  float f4; 
	from_vax_r4( &f, &f4, &n );
  
	return(f4);
}

double vax2d (double le_dbl)
{
  int n =1;
	double d8;
	
	from_vax_d8( (double *)&le_dbl, (double *)&d8, (int *)&n );	 
  return(d8);
}

/* ------------------------------------------------------------------------------ */


int is_big_endian(void)
{
    
    int iFlagBigEndian = -1;
    union {
        short s;
        char c[sizeof(short)];
    } un;
/*    un.s = Endian_Word_Conversion((0x0102));*/
    un.s = 0x0102;
    if(sizeof(short) == 2)
    {
        if(un.c[0] == 1 && un.c[1] == 2){
            iFlagBigEndian = 1; 
        }
        else if(un.c[0] == 2 && un.c[1] == 1){
            iFlagBigEndian = 0; 
        }
    }
    return(iFlagBigEndian);
}


/* General byte swapping routine */
 void SwapBytes(void *source, int size)
{
    typedef unsigned char TwoBytes[2];
    typedef unsigned char FourBytes[4];
    typedef unsigned char EightBytes[8];

    unsigned char temp;

    if(size == 2)
    {
        TwoBytes *src = (TwoBytes *)source;
        temp = (*src)[0];
        (*src)[0] = (*src)[1];
        (*src)[1] = temp;

        return;
    }

    if(size == 4)
    {
        FourBytes *src = (FourBytes *)source;
        temp = (*src)[0];
        (*src)[0] = (*src)[3];
        (*src)[3] = temp;

        temp = (*src)[1];
        (*src)[1] = (*src)[2];
        (*src)[2] = temp;

        return;
    }

    if(size == 8)
    {
        EightBytes *src = (EightBytes *)source;
        temp = (*src)[0];
        (*src)[0] = (*src)[7];
        (*src)[7] = temp;

        temp = (*src)[1];
        (*src)[1] = (*src)[6];
        (*src)[6] = temp;

        temp = (*src)[2];
        (*src)[2] = (*src)[5];
        (*src)[5] = temp;

        temp = (*src)[3];
        (*src)[3] = (*src)[4];
        (*src)[4] = temp;

        return;
    }

}
