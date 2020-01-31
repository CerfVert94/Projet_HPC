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
#include "mynrutil.h"

void memset_ui8matrix(uint8 **X, uint8 value, long nrl, long nrh, long ncl, long nch)
/* ------------------------------------------------------------------------- */
{ 
      memset((uint8*)(&X[nrl][ncl]- NR_END), value, sizeof(**X) * (nrh - nrl + 1) * (nch - ncl + 1));
}

int memcmp_ui8matrix(uint8 **X, uint8 **Y, long nrl, long nrh, long ncl, long nch)
/* ------------------------------------------------------------------------- */
{ 
    return memcmp(X[nrl] + ncl, Y[nrl] + ncl, sizeof(uint8) * (nrh - nrl + 1) * (nch - ncl + 1));
}

