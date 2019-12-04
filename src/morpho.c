#include "nrdef.h"
#include "nrutil.h"
#include <img.h>
#include <stdlib.h>
#include <stdio.h>
//#include <malloc.h>
#include <morpho.h>
#include <test_morpho.h>
#include <util.h>

extern const char * nom_func;
p_struct_elem_dim compute_struct_elem_dim(long x0, long y0, long nrow, long ncol)
{
	p_struct_elem_dim s;

	s = (p_struct_elem_dim) malloc(sizeof(struct_elem_dim));
	
	if (!s) 
		exit_on_error("malloc failed");
		
	// Define the dimension.
	s->nrow = nrow; s->ncol = ncol;
	// Define the origin.
	s->x0= x0;  s->y0= y0;
	// Compute the lower / upper limits
	s->nrl =  y0 - (nrow - 1); 
	s->nrh = -y0 + (nrow - 1);
	s->ncl =  x0 - (ncol - 1); 
	s->nch = -x0 + (ncol - 1);
	return s;
}
void free_structuring_element(p_struct_elem_dim s)
{
	free(s);
}

inline void ui8matrix_sequence_naive(uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
	uint8 **ppPreOutput0, **ppPreOutput1;
	ppPreOutput0 = ui8matrix(nrl + s->nrl, nrh + s->nrh, ncl + s->ncl, nch + s->nch);
	ppPreOutput1 = ui8matrix(nrl + s->nrl, nrh + s->nrh, ncl + s->ncl, nch + s->nch);
	ui8matrix_erosion_naive (ppInput     , nrl, nrh, ncl, nch, s, ppPreOutput0);
	ui8matrix_dilation_naive(ppPreOutput0, nrl, nrh, ncl, nch, s, ppPreOutput1);
	ui8matrix_dilation_naive(ppPreOutput1, nrl, nrh, ncl, nch, s, ppPreOutput0);
	ui8matrix_erosion_naive (ppPreOutput0, nrl, nrh, ncl, nch, s, ppOutput);
	free_ui8matrix(ppPreOutput0, nrl + s->nrl, nrh + s->nrh, ncl + s->ncl, nch + s->nch);
	free_ui8matrix(ppPreOutput1, nrl + s->nrl, nrh + s->nrh, ncl + s->ncl, nch + s->nch);
}
void ui8matrix_erosion_naive(uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
	long row, col, x, y;
	long snrl = s->nrl, snrh = s->nrh, sncl = s->ncl, snch = s->nch;

	// Erode
	for (row = nrl; row < nrh + 1; row++){
		for (col = ncl; col < nch + 1; col++) {
			// Copy the input to the output
			ppOutput[row][col] = ppInput[row][col];

			// Apply morpho to the output
			for (y = snrl; y < snrh + 1; y++ ) 
				for (x = sncl; x < snch + 1; x++ ) 
					ppOutput[row][col] &= ppInput[row + y][col + x];
		}
	}
}
void ui8matrix_dilation_naive(uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
	long row, col, x, y;
	long snrl = s->nrl, snrh = s->nrh, sncl = s->ncl, snch = s->nch;
	// dilate
	for (row = nrl; row < nrh + 1; row++){
		for (col = ncl; col < nch + 1; col++) {
			// Copy the input to the output
			ppOutput[row][col] = ppInput[row][col];
			
			// Apply morpho to the output
			for (y = snrl; y < snrh + 1; y++ ) 
				for (x = sncl; x < snch + 1; x++ ) 
					ppOutput[row][col] |= ppInput[row + y][col + x];
		}
	}
}



#define _IMG_DIMENSION(img, s) img->nrl - s->nrl, img->nrh - s->nrh, img->ncl - s->ncl, img->nch - s->nch
// #define _IMG_DIMENSION(img) img->nrl, img->nrh, img->ncl, img->nch
void image_chain_processing(p_image img, p_struct_elem_dim s, int idx)
{
	uint8 **ppOutput;
	char filename[128];
	long nrl, nrh, ncl, nch;

	ppOutput = ui8matrix(img->nrl, img->nrh, img->ncl, img->nch);
	nrl = img->nrl + BORD;//s->nrl;
	nrh = img->nrh - BORD;//s->nrh;
	ncl = img->ncl + BORD;//s->ncl;
	nch = img->nch - BORD;//s->nch;

	ui8matrix_erosion_naive (img->E, nrl, nrh, ncl, nch, s, ppOutput);	
	ui8matrix_dilation_naive(ppOutput,  nrl, nrh, ncl, nch, s, img->Omega);	
	ui8matrix_dilation_naive(img->Omega, nrl, nrh, ncl, nch, s, ppOutput);	
	ui8matrix_erosion_naive (ppOutput, nrl, nrh, ncl, nch, s, img->Omega);	

	for (int i = img->nrl; i < img->nrh; i++) {
		for (int j = img->ncl; j < img->nch; j++) {
			// printf("%d", img->E[i][j]);
			img->Omega[i][j] *= 255;
		}
		// putchar('\n');
	}
	// // Save the image (debug)
	sprintf(filename, "../car3_sigma/car_%d.pgm", idx);
	SavePGM_ui8matrix(img->Omega, nrl, nrh, ncl, nch, filename);
	free_ui8matrix(ppOutput, img->nrl, img->nrh, img->ncl, img->nch);
}


// Without optimisation
uint8 dilation_naive(uint8** ppInput, long row, long col, p_struct_elem_dim s)
{
	uint8 pixel = ppInput[row][col];
	// Dilate
	for (long y = s->nrl; y < s->nrh + 1; y++ ) {
		for (long x = s->ncl; x < s->nch + 1; x++ ) {
			pixel |= ppInput[row + y][col + x];
		}
	}
	return pixel;
}
uint8 erosion_naive(uint8** ppInput, long row, long col, p_struct_elem_dim s)
{
	uint8 pixel = ppInput[row][col];
	// Erode
	for (long y = s->nrl; y < s->nrh + 1; y++ ) {
		for (long x = s->ncl; x < s->nch + 1; x++ ) {
			pixel &= ppInput[row + y][col + x];
		}
	}
	return pixel;
}
