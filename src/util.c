#include <nrdef.h>
#include <stdio.h>
#include <stdlib.h>
#include <util.h>
void binary_to_octal_ui8matrix(uint8 **ppInput, long nrl, long nrh, long ncl, long nch)
{
    for (long i = nrl; i < nrh + 1; i++){
        for (long j = ncl; j < nch + 1; j++){
            if (ppInput[i][j] > 0) 
                ppInput[i][j] = 255;
            else 
                ppInput[i][j] = 0;
        }
    }
}
void octal_to_binary_ui8matrix(uint8 **ppInput, long nrl, long nrh, long ncl, long nch)
{
    for (long i = nrl; i < nrh + 1; i++){
        for (long j = ncl; j < nch + 1; j++){
            if (ppInput[i][j] > 0) 
                ppInput[i][j] = 1;
            else 
                ppInput[i][j] = 0;
        }
    }
}