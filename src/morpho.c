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
	if (!s) {
		perror("malloc failed");
		exit(EXIT_FAILURE);
	}
	// Define dimension.
	s->nrow = nrow; s->ncol = ncol;
	// Define origin.
	s->x0= x0;  s->y0= y0;
	// Compute the size of borders
	s->nrl = -y0; s->nrh = (nrow - 1) - y0;
	s->ncl = -x0; s->nch = (ncol - 1) - x0;
	return s;
}
void free_structuring_element(p_struct_elem_dim s)
{
	free(s);
}

void ui8matrix_erosion_naive(uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
	long row, col, x, y;
	// Erode
	for (row = nrl; row < nrh + 1; row++){
		for (col = ncl; col < nch + 1; col++) {
			// Copy the input to the output
			ppOutput[row][col] = ppInput[row][col];

			// Apply morpho to the output
			for (y = s->nrl; y < s->nrh + 1; y++ ) 
				for (x = s->ncl; x < s->nch + 1; x++ ) 
					ppOutput[row][col] &= ppInput[row + y][col + x];
		}
	}
}
void ui8matrix_dilation_naive(uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput)
{
	long row, col, x, y;
	// dilate
	for (row = nrl; row < nrh + 1; row++){
		for (col = ncl; col < nch + 1; col++) {
			// Copy the input to the output
			ppOutput[row][col] = ppInput[row][col];
			
			// Apply morpho to the output
			for (y = s->nrl; y < s->nrh + 1; y++ ) 
				for (x = s->ncl; x < s->nch + 1; x++ ) 
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


/*
Test for morpho
	A. Erosion (3x3 / 5x5)
		a. Produce 1
			i. Full rectangle 
		b. Produce 0
			i.  Empty rectangle
			ii. Test corners
			iii. Test edges
	B. Dilation
		a. Produce 1
			i. Full rectangle 
			ii. Test corners
			iii. Test edges
		b. Produce 0
			i.  Empty rectangle
*/

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
