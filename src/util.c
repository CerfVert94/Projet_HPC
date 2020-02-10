#include <nrdef.h>
#include <nrutil.h>
#include <mynrutil.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <util.h>
#include <assert.h>

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
            if (ppInput[i][j] > 127) 
                ppInput[i][j] = 1;
            else 
                ppInput[i][j] = 0;
        }
    }
}


uint8 **vpack_binary_ui8matrix  (uint8 **X, long nrl, long nrh, long ncl, long nch, long *packed_nrl, long *packed_nrh, long *packed_ncl, long *packed_nch)
{
    long packed_nrow, packed_ncol;
    uint8** Y;
    long row, col, row_prime, cnt;

    (*packed_nrl) = pack8(nrl);
    (*packed_nrh) = pack8(nrh);
    (*packed_ncl) = ncl;// < 0 ? ncl / 8 + ROUNDUP_OVER8(ncl % 8) : (ncl + 1) / 8 + ROUNDUP_OVER8((ncl + 1) % 8);
    (*packed_nch) = nch;// < 0 ? nch / 8 + ROUNDUP_OVER8(nch % 8) : (nch + 1) / 8 + ROUNDUP_OVER8((nch + 1) % 8);

    Y = ui8matrix(*packed_nrl, *packed_nrh, *packed_ncl, *packed_nch);
    memset_ui8matrix(Y, 0, *packed_nrl, *packed_nrh, *packed_ncl, *packed_nch);
    
 
    return Y;
}
uint8 **unvpack_binary_ui8matrix(uint8 **X, long nrl, long nrh, long ncl, long nch)
{
    long packed_nrow, packed_ncol;
    uint8** Y;
    long row, col, row_prime, cnt;
    long packed_nrl, packed_nrh, packed_ncl, packed_nch;

    packed_nrl = pack8(nrl);
    packed_nrh = pack8(nrh);
    packed_ncl = ncl;// < 0 ? ncl / 8 + ROUNDUP_OVER8(ncl % 8) : (ncl + 1) / 8 + ROUNDUP_OVER8((ncl + 1) % 8);
    packed_nch = nch;// < 0 ? nch / 8 + ROUNDUP_OVER8(nch % 8) : (nch + 1) / 8 + ROUNDUP_OVER8((nch + 1) % 8);

    Y = ui8matrix(nrl, nrh, ncl, nch);
    memset_ui8matrix(Y, 0, nrl, nrh, ncl, nch);
    for (col = ncl; col < nch + 1; col++) {
        row_prime = packed_nrl;
        for (row = nrl; row < nrh + 1; row++) {
            Y[row][col] = (X[row_prime][col] >> cnt) & 0x1;
            cnt = (cnt + 1) % 8;
            if (cnt == 0)  ++row_prime;
        } 
    } 
    return Y;
}

#define SE_NRL -1
#define SE_NRH  1
#define SE_NCL -1
#define SE_NCH  1

void packing_test(char *filename, pack_func_t pack_func, unpack_func_t unpack_func, const char *func_name, bool display)
{
      
    long nrl, ncl, nrh, nch, row, col;
    long packed_nrl, packed_ncl, packed_nrh, packed_nch, bord;
    uint8 **image, **X, **Y, **Z;

    image = LoadPGM_ui8matrix(filename, &nrl, &nrh, &ncl, &nch);
    octal_to_binary_ui8matrix(image, nrl, nrh, ncl, nch);
    // Test Input / Test Output
    X = ui8matrix(nrl, nrh, ncl, nch);
    // Full zero intialization & copy image 
    memset_ui8matrix(X, 0, nrl, nrh, ncl, nch);
    copy_ui8matrix_ui8matrix(image, nrl, nrh, ncl, nch, X);

    printf("Pack / Unpack test : "LALIGNED_STR" (%-30s)\n", func_name, filename);
    printf("\tTesting for %ld x %ld\n", nch + 1, nrh + 1);
    // Packed Output
    Y = fcpacked_ui8matrix(nrl, nrh, ncl, nch, &packed_nrl, &packed_nrh, &packed_ncl, &packed_nch, &bord);

    // pack_func(X, 
    // Unpacked Output
    // Z = unpack_func(Y, nrl, nrh, ncl, nch, packed_nrh, packed_ncl, packed_nch, bord);

    if (display) {
        display_ui8matrix(X, nrl, nrh, ncl, nch, "%u", "Valid output");
        display_ui8matrix(Y, packed_nrl, packed_nrh, packed_ncl, packed_nch, "%u", "Packed");
        display_ui8matrix(Z, nrl, nrh, ncl, nch, "%u", "Unpacked");
    }
    assert(!memcmp(X[nrl] + ncl, Z[nrl] + ncl, (nrh - nrl + 1) * (nch - ncl + 1) * sizeof(**Z)));
    puts("Passed the test.\n");

    free_ui8matrix(Y, packed_nrl, packed_nrh, packed_ncl, packed_nch); 
    free_ui8matrix(Z, nrl, nrh, ncl, nch); 
    free_ui8matrix(X, nrl, nrh, ncl, nch);
	free_ui8matrix(image, nrl, nrh, ncl, nch);
}

#define BORDER_SIZE 2
uint8 **fcpacked_ui8matrix(long nrl, long nrh, long ncl, long nch, long *packed_nrl, long *packed_nrh, long *packed_ncl, long *packed_nch, long *bord)
{
    uint8** X;
    *packed_nrl =       nrl ;
    *packed_nrh =       nrh ;
    *packed_ncl = pack8(ncl);
    *packed_nch = pack8(nch);
    *bord       = BORDER_SIZE;
    X = ui8matrix(*packed_nrl - *bord, *packed_nrh + *bord, *packed_ncl - *bord, *packed_nch + *bord);
    memset_ui8matrix(X, 0, *packed_nrl - *bord, *packed_nrh + *bord, *packed_ncl - *bord, *packed_nch + *bord);
    return X;
}
uint8 **frpacked_ui8matrix(long nrl, long nrh, long ncl, long nch, long *packed_nrl, long *packed_nrh, long *packed_ncl, long *packed_nch, long *bord)
{
    uint8** X;

    *packed_nrl = pack8(nrl);
    *packed_nrh = pack8(nrh);
    *packed_ncl =       ncl ;
    *packed_nch =       nch ;
    *bord       = BORDER_SIZE;
    X = ui8matrix(*packed_nrl - *bord, *packed_nrh + *bord, *packed_ncl - *bord, *packed_nch + *bord);
    memset_ui8matrix(X, 0, *packed_nrl - *bord, *packed_nrh + *bord, *packed_ncl - *bord, *packed_nch + *bord);
    return X;
}
uint8 **hcpacked_ui8matrix(long nrl, long nrh, long ncl, long nch, long *packed_nrl, long *packed_nrh, long *packed_ncl, long *packed_nch, long *bord)
{
    uint8** X;

    *packed_nrl = pack8(nrl);
    *packed_nrh = pack8(nrh);
    *packed_ncl =          ncl ;
    *packed_nch =          nch ;
    *bord       = BORDER_SIZE;
    X = ui8matrix(*packed_nrl - *bord, *packed_nrh + *bord, *packed_ncl - *bord, *packed_nch + *bord);
    memset_ui8matrix(X, 0, *packed_nrl - *bord, *packed_nrh + *bord, *packed_ncl - *bord, *packed_nch + *bord);
    return X;
}
uint8 **hrpacked_ui8matrix(long nrl, long nrh, long ncl, long nch, long *packed_nrl, long *packed_nrh, long *packed_ncl, long *packed_nch, long *bord)
{
    uint8** X;

    *packed_nrl = pack8(nrl);
    *packed_nrh = pack8(nrh);
    *packed_ncl =          ncl ;
    *packed_nch =          nch ;
    *bord       = BORDER_SIZE;
    X = ui8matrix(*packed_nrl - *bord, *packed_nrh + *bord, *packed_ncl - *bord, *packed_nch + *bord);
    memset_ui8matrix(X, 0, *packed_nrl - *bord, *packed_nrh + *bord, *packed_ncl - *bord, *packed_nch + *bord);
    return X;
}

void free_packed_ui8matrix (uint8 **X, long  packed_nrl, long  packed_nrh, long  packed_ncl, long  packed_nch, long  bord)
{
    free_ui8matrix(X, packed_nrl - bord, packed_nrh + bord, packed_ncl - bord, packed_nch + bord);
}
void fcpack_ui8matrix_ui8matrix (uint8 **X, long nrl, long nrh, long ncl, long nch, long packed_nrl, long packed_nrh, long packed_ncl, long packed_nch, long bord, uint8 **Y)
{
    long row, col, col_prime, ncl_prime, nch_prime, cnt;
    uint8 init_val, final_val;

    packed_nrl -= bord;
    packed_nrh += bord;
    packed_ncl -= bord;
    packed_nch += bord;

    ncl_prime = packed_ncl * 8;
    nch_prime = packed_nch * 8;

    
    for (row = nrl; row < nrh + 1; row++) {
        col_prime = packed_ncl;
        cnt = 0;
        Y[row][col_prime] = 0;

        init_val  = 0;//X[row][ncl];
        final_val = 0;//X[row][nch];
        for (col = ncl_prime; col < nch_prime + 1; col++) {
            
                 if (ncl <= col && col <= nch) Y[row][col_prime] |= X[row][col] << (7 - (cnt % 8));
            else if (col < ncl)                Y[row][col_prime] |= init_val    << (7 - (cnt % 8));
            else if (col > nch)                Y[row][col_prime] |= final_val   << (7 - (cnt % 8));

            cnt++;
            if (cnt % 8 == 0 && 0 < cnt)       Y[row][++col_prime] = 0;
        }
    }
    // Prologue
    for (row = packed_nrl; row < nrl + 1; row++) {

    }
    // Epilouge
    for (row = nrh; row < packed_nch + 1; row++) {

    }
}

void unfcpack_ui8matrix_ui8matrix (uint8 **X, long nrl, long nrh, long ncl, long nch, long packed_nrl, long packed_nrh, long packed_ncl, long packed_nch, long bord, uint8 **Y)
{
    long row, col, col_prime,  ncl_prime, nch_prime, cnt;
     
    packed_ncl -= bord;
    packed_nch += bord;

    ncl_prime = packed_ncl * 8;
    nch_prime = packed_nch * 8;

    for (row = nrl; row < nrh + 1; row++) {
        col_prime = packed_ncl;
        cnt = 0;

        for (col = ncl_prime; col < nch_prime + 1; col++) {
            if (ncl <= col && col <= nch) Y[row][col] = (X[row][col_prime] >> (7 - (cnt % 8))) & 0x1;

            cnt++;
            if (cnt % 8 == 0 && 0 < cnt) col_prime++;
        }
    }
}
