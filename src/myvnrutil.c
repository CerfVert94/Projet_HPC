
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <xmmintrin.h>
#include <emmintrin.h>
#include <mmintrin.h>

#include "nrdef.h"
#include "vnrdef.h"
#include "vnrutil.h"
#include "myvnrutil.h"
#include "mynrutil.h"

uint8** vui8matrix_to_ui8matrix(vuint8** vX, long i0, long i1, int j0, int j1, long *nrl, long *nrh, long *ncl, long *nch)
/*--------------------------------------------------------------------------------*/
{
    int i, j, card;
    vuint8 x;

    card = card_vuint8();
    *nrl = i0;        *nrh = i1;
    *ncl = j0 * card; *nch = (j0 + (j1 - j0 + 1)) * card - 1;
    
    uint8 **Y = filled_ui8matrix(*nrl, *nrh, *ncl, *nch, 0);

    
    for(i=i0; i<=i1; i++) {
        for(j=j0; j<=j1; j++) {
            printf("j : %ld\n",j);
            vuint8 T[1];
            uint8 *U;
            x = _mm_load_si128(&vX[i][j]);
            _mm_store_si128(T, x);                  
            memmove(&Y[i][j * card], (uint8*)T, sizeof(uint8) * card);

        }
    }
    return Y;
}