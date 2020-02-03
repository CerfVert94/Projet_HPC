#include "nrdef.h"
#include "nrutil.h"
#include <img.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
//#include <malloc.h>
#include "vnrdef.h"
#include "vnrutil.h"
#include <img_SIMD.h>
#include <util.h>
#include <morpho.h>


#define SE_NRL -1
#define SE_NRH  1
#define SE_NCL -1
#define SE_NCH  1


extern const char * nom_func;

inline void ui8matrix_sequence_naive(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	// uint8 **ppPreOutput0, **ppPreOutput1;
	// ppPreOutput0 = ui8matrix(nrl + SE_NRL, nrh + SE_NRH, ncl + SE_NCL, nch + SE_NCH);
	// ppPreOutput1 = ui8matrix(nrl + SE_NRL, nrh + SE_NRH, ncl + SE_NCL, nch + SE_NCH);
	// memset_ui8matrix(ppPreOutput0, 0, nrl + SE_NRL, nrh + SE_NRH, ncl + SE_NCL, nch + SE_NCH);
	// ui8matrix_erosion_naive (X     , nrl, nrh, ncl, nch, ppPreOutput0);
	// memset_ui8matrix(ppPreOutput1, 0, nrl + SE_NRL, nrh + SE_NRH, ncl + SE_NCL, nch + SE_NCH);
	// ui8matrix_dilation_naive(ppPreOutput0, nrl, nrh, ncl, nch, ppPreOutput1);
	// memset_ui8matrix(ppPreOutput0, 0, nrl + SE_NRL, nrh + SE_NRH, ncl + SE_NCL, nch + SE_NCH);
	// ui8matrix_dilation_naive(ppPreOutput1, nrl, nrh, ncl, nch, ppPreOutput0);
	// ui8matrix_erosion_naive (ppPreOutput0, nrl, nrh, ncl, nch, Y);
	// free_ui8matrix(ppPreOutput0, nrl + SE_NRL, nrh + SE_NRH, ncl + SE_NCL, nch + SE_NCH);
	// free_ui8matrix(ppPreOutput1, nrl + SE_NRL, nrh + SE_NRH, ncl + SE_NCL, nch + SE_NCH);
}
void ui8matrix_erosion_naive(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	long row, col, x, y;
	long snrl = SE_NRL, snrh = SE_NRH, sncl = SE_NCL, snch = SE_NCH;

	// Erode
	for (row = nrl; row < nrh + 1; row++){
		for (col = ncl; col < nch + 1; col++) {
			// Copy the input to the output
			Y[row][col] = X[row][col];

			// Apply morpho to the output
			for (y = snrl; y < snrh + 1; y++ ) 
				for (x = sncl; x < snch + 1; x++ )  
					Y[row][col] &= X[row + y][col + x];
		}
	}
}
void ui8matrix_dilation_naive(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	long row, col, x, y;
	long snrl = SE_NRL, snrh = SE_NRH, sncl = SE_NCL, snch = SE_NCH;
	// dilate
	for (row = nrl; row < nrh + 1; row++){
		for (col = ncl; col < nch + 1; col++) {
			// Copy the input to the output
			Y[row][col] = X[row][col];
			
			// Apply morpho to the output
			for (y = snrl; y < snrh + 1; y++ )
				for (x = sncl; x < snch + 1; x++ ) 
					Y[row][col] |= X[row + y][col + x];
		}
	}
}



#define _IMG_DIMENSION(img, s) img->nrl - SE_NRL, img->nrh - SE_NRH, img->ncl - SE_NCL, img->nch - SE_NCH
// #define _IMG_DIMENSION(img) img->nrl, img->nrh, img->ncl, img->nch
void image_chain_processing(p_image img, int idx)
{
	// uint8 **Y;
	// char filename[128];
	// long nrl, nrh, ncl, nch;

	// Y = ui8matrix(img->nrl, img->nrh, img->ncl, img->nch);
	// nrl = img->nrl + BORD;//SE_NRL;
	// nrh = img->nrh - BORD;//SE_NRH;
	// ncl = img->ncl + BORD;//SE_NCL;
	// nch = img->nch - BORD;//SE_NCH;

	// ui8matrix_erosion_naive (img->E, nrl, nrh, ncl, nch, Y);	
	// ui8matrix_dilation_naive(Y,  nrl, nrh, ncl, nch, img->Omega);	
	// ui8matrix_dilation_naive(img->Omega, nrl, nrh, ncl, nch, Y);	
	// ui8matrix_erosion_naive (Y, nrl, nrh, ncl, nch, img->Omega);	

	// for (int i = img->nrl; i < img->nrh; i++) {
	// 	for (int j = img->ncl; j < img->nch; j++) {
	// 		// printf("%d", img->E[i][j]);
	// 		img->Omega[i][j] *= 255;
	// 	}
	// 	// putchar('\n');
	// }
	// // // Save the image (debug)
	// sprintf(filename, "../car3_sigma/car_%d.pgm", idx);
	// SavePGM_ui8matrix(img->Omega, nrl, nrh, ncl, nch, filename);
	// free_ui8matrix(Y, img->nrl, img->nrh, img->ncl, img->nch);
}


// Without optimisation
uint8 dilation_naive(uint8** X, long row, long col)
{
	uint8 pixel = X[row][col];
	// Dilate
	for (long y = SE_NRL; y < SE_NRH + 1; y++ ) {
		for (long x = SE_NCL; x < SE_NCH + 1; x++ ) {
			pixel |= X[row + y][col + x];
		}
	}
	return pixel;
}
uint8 erosion_naive(uint8** X, long row, long col)
{
	uint8 pixel = X[row][col];
	// Erode
	for (long y = SE_NRL; y < SE_NRH + 1; y++ ) {
		for (long x = SE_NCL; x < SE_NCH + 1; x++ ) {
			pixel &= X[row + y][col + x];
		}
	}
	return pixel;
}
