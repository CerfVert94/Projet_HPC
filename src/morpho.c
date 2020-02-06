#include "nrdef.h"
#include "nrutil.h"
#include "mynrutil.h"
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
	memset_ui8matrix(temp_buffer, 0, nrl-2, nrh+2, ncl-2, nch+2);
	
	
	ui8matrix_erosion_naive  (X, nrl, nrh, ncl, nch, NULL, temp_buffer); memset_ui8matrix(Y, 0, nrl-2, nrh+2, ncl-2, nch+2);
	ui8matrix_dilation_naive (temp_buffer, nrl, nrh, ncl, nch, NULL, Y); memset_ui8matrix(temp_buffer, 0, nrl-2, nrh+2, ncl-2, nch+2);
	ui8matrix_dilation_naive (Y, nrl, nrh, ncl, nch, NULL, temp_buffer); memset_ui8matrix(Y, 0, nrl-2, nrh+2, ncl-2, nch+2);
	ui8matrix_erosion_naive  (temp_buffer, nrl, nrh, ncl, nch, NULL, Y);
	// memcpy_ui8matrix(temp_buffer, nrl - 2, nrh + 2, ncl - 2, nrh + 2, Y);
	
	// display_ui8matrix(Y, nrl-2, nrh+2, ncl-2, nch+2, "%4u", "Sequence");
	
	
	
	// display_ui8matrix(Y          , nrl-2, nrh+2, ncl-2, nch+2, "%4u", "3rd Erosion (Naive)");
}
void ui8matrix_erosion_naive(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	long row, col, x, y;
	long snrl = -1, snrh = 1, sncl = -1, snch = 1;

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
	long snrl = -1, snrh = 1, sncl = -1, snch = 1;
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

void ui8matrix_erosion5_naive(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	long row, col, x, y;
	long snrl = -2, snrh = 2, sncl = -2, snch = 2;

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
void ui8matrix_dilation5_naive(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	long row, col, x, y;
	long snrl = -2, snrh = 2, sncl = -2, snch = 2;
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
