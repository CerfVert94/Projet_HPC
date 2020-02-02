/* -------------- */
/* --- nrutil --- */
/* -------------- */

/*
 * Copyright (c) 2000 - 2007, Lionel Lacassagne
 * Ensta version
 */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdarg.h>
#include <ctype.h> 
#include <string.h>
#include <math.h> /* fabs */

#include "nrdef.h"
#include "nrutil.h"
#include "mynrutil.h"
/* ------------------------------------------------------------------------- */
void memset_ui8matrix(uint8 **X, uint8 value, long nrl, long nrh, long ncl, long nch)
/* ------------------------------------------------------------------------- */
{ 
      memset((uint8*)(&X[nrl][ncl]- NR_END), value, sizeof(**X) * (nrh - nrl + 1) * (nch - ncl + 1));
}
/* ------------------------------------------------------------------------- */
void memcpy_ui8matrix(uint8 **X, long nrl, long nrh, long ncl, long nch, uint8 **Y)
/* ------------------------------------------------------------------------- */
{ 
      memcpy((uint8*)(&Y[nrl][ncl]- NR_END),(uint8*)(&X[nrl][ncl]- NR_END),  sizeof(uint8) * (nrh - nrl + 1) * (nch - ncl + 1));
}
/* ------------------------------------------------------------------------- */
int memcmp_ui8matrix(uint8 **X, uint8 **Y, long nrl, long nrh, long ncl, long nch)
/* ------------------------------------------------------------------------- */
{ 
    return memcmp(X[nrl] + ncl, Y[nrl] + ncl, sizeof(uint8) * (nrh - nrl + 1) * (nch - ncl + 1));
}
/* ------------------------------------------------------------------------- */
uint8** filled_ui8matrix(long nrl, long nrh, long ncl, long nch, uint8 value)
/* ------------------------------------------------------------------------- */
{
    uint8 **X = ui8matrix(nrl, nrh, ncl, nch);
    memset_ui8matrix(X, value, nrl, nrh, ncl, nch);
    return X;
}
