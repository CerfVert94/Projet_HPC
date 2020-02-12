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

void ui8matrix_sequence_naive(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **Y, uint8 **Z)
{
	long row, col, x, y;
	long snrl = -1, snrh = 1, sncl = -1, snch = 1;

	memset_ui8matrix(Y, 0, nrl-2, nrh+2, ncl-2, nch+2);
	memset_ui8matrix(Z, 0, nrl-2, nrh+2, ncl-2, nch+2);
	// for (long row = nrl; row < nrh + 1; row++){
	// 	X[row][ncl - 1] = X[row][ncl];
	// 	X[row][nch + 1] = X[row][nch];
	// }
	// display_ui8matrix(X, nrl - 2, nrh + 2, ncl - 2, nch + 2, "%4u", "INPUT");
	ui8matrix_erosion_naive  (X, nrl, nrh, ncl, nch, NULL, Y); memset_ui8matrix(Z, 0, nrl-2, nrh+2, ncl-2, nch+2);
	// display_ui8matrix(Y, nrl - 2, nrh + 2, ncl - 2, nch + 2, "%4u", "E3 (Naive)");
	// display_ui8matrix(Y, nrl - 2, nrh + 2, -2 + ncl, 32 + 2, "%4u", "Erosion seq_naive");


	// for (long row = nrl; row < nrh + 1; row++){
	// 	Y[row][ncl - 1] = Y[row][ncl];
	// 	Y[row][nch + 1] = Y[row][nch];
	// }
	ui8matrix_dilation_naive (Y, nrl, nrh, ncl, nch, NULL, Z); memset_ui8matrix(Y, 0, nrl-2, nrh+2, ncl-2, nch+2);
	// for (long row = nrl; row < nrh + 1; row++){
	// 	Z[row][ncl - 1] = Z[row][ncl];
	// 	Z[row][nch + 1] = Z[row][nch];
	// }
	ui8matrix_dilation_naive (Z, nrl, nrh, ncl, nch, NULL, Y); memset_ui8matrix(Z, 0, nrl-2, nrh+2, ncl-2, nch+2);														  
	// display_ui8matrix(Y, nrl - 2, nrh + 2, ncl - 2, nch + 2, "%4u", "E3 - D5 (Naive)");
	// for (long row = nrl; row < nrh + 1; row++){
	// 	Y[row][ncl - 1] = Y[row][ncl];
	// 	Y[row][nch + 1] = Y[row][nch];
	// }
	// display_ui8matrix(Y, nrl - 2, nrh + 2, ncl - 2, nch + 2, "%4u", "E3-D5 (Naive)");
	// display_ui8matrix(Y, nrl - 2, nrh + 2, -2 + ncl, nch + 2, "%4u", "E3-D5 (Naive)");
	ui8matrix_erosion_naive(Y, nrl, nrh, ncl, nch, NULL, Z);
	// display_ui8matrix(Z, nrl - 2, nrh + 2, ncl - 2, nch + 2, "%4u", "E3-D5-E3 (Naive)");
	
	
	
	display_ui8matrix(Z, nrl - 2, nrh + 2, ncl, nch , "%4u", "EDDE");
	// getchar();
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
