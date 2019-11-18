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
p_struct_elem_dim compute_struct_elem_dim(long orix, long oriy, long nrow, long ncol)
{
	p_struct_elem_dim s;

	s = (p_struct_elem_dim) malloc(sizeof(struct_elem_dim));
	// Define dimension.
	s->nrow = nrow; s->ncol = ncol;
	// Define origin.
	s->orix= orix;  s->oriy= oriy;
	// Compute the size of borders
	s->nthbord = oriy; s->nbhbord = (nrow - 1) - oriy;
	s->nlvbord = orix; s->nrvbord = (ncol - 1) - orix;
}
void free_structuring_element(p_struct_elem_dim s)
{
	free(s);
}

// Without optimisation
uint8** ui8matrix_dilation(uint8** input, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s)
{
	uint8 **bord_input, **output;
	long col, row, y, x;
	uint8 pixel;
	
	// Create a new matrix with edges / output matrix
	bord_input  = ui8matrix(nrl - s->nthbord, nrh + s->nbhbord, ncl - s->nlvbord, nch + s->nrvbord);
	output 		= ui8matrix(nrl, nrh, ncl, nch);
	
	// Copy the original image to the input matrix.
	copy_ui8matrix_ui8matrix(input, nrl, nrh, ncl, nch, bord_input);
	
	// Erode
	for (row = nrl; row < nrh + 1; row++){
		for (col = ncl; col < nch + 1; col++) {
			pixel = bord_input[row][col];
			for (y = -s->nthbord; y < s->nbhbord + 1; y++ ) {
				for (x = -s->nlvbord; x < s->nrvbord + 1; x++ ) {
					pixel |= bord_input[row + y][col + x];
				}
			}
			output[row][col] = pixel;
		}
	}
	
	// Free the bordered matrix.
	free_ui8matrix(bord_input, nrl - s->nthbord, nrh + s->nbhbord, ncl - s->nlvbord, nch + s->nrvbord);
	return output;
}
uint8** ui8matrix_erosion(uint8** input, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s)
{
	uint8 **bord_input, **output;
	long col, row, y, x;
	uint8 pixel;
	
	// Create a new matrix with edges / output matrix
	bord_input  = ui8matrix(nrl - s->nthbord, nrh + s->nbhbord, ncl - s->nlvbord, nch + s->nrvbord);
	output 		= ui8matrix(nrl, nrh, ncl, nch);

	// Copy the original image to the input matrix.
	copy_ui8matrix_ui8matrix(input, nrl, nrh, ncl, nch, bord_input);	
	// Dilate
	for (row = nrl; row < nrh + 1; row++){
		for (col = ncl; col < nch + 1; col++) {
			pixel = bord_input[row][col];
			for (y = -s->nthbord; y < s->nbhbord + 1; y++ ) {
				for (x = -s->nlvbord; x < s->nrvbord + 1; x++ ) {
					pixel &= bord_input[row + y][col + x];
				}
			}
			output[row][col] = pixel;
		}
	}
	// Free the bordered matrix.
	free_ui8matrix(bord_input, nrl - s->nthbord, nrh + s->nbhbord, ncl - s->nlvbord, nch + s->nrvbord);
	return output;
}


void image_dilation(p_image img, p_struct_elem_dim s)
{
	uint8 **output;
	long nrl, nrh, ncl, nch;
	
	nrl = img->nrl; nrh = img->nrh;
	ncl = img->ncl;	nch = img->nch;

	output = ui8matrix_dilation(img->I, nrl, nrh, ncl, nch, s);	
	copy_ui8matrix_ui8matrix(output, nrl, nrh, ncl, nch, img->Omega);
	// Save the image (debug)
	SavePGM_ui8matrix(img->Omega, nrl, nrh, ncl, nch,"output_dilation.pgm");
	free_ui8matrix(output, nrl, nrh, ncl, nch);
}
void image_erosion(p_image img, p_struct_elem_dim s)
{
	uint8 **output;
	long nrl, nrh, ncl, nch;
	
	nrl = img->nrl; nrh = img->nrh;
	ncl = img->ncl;	nch = img->nch;

	output = ui8matrix_erosion(img->I, nrl, nrh, ncl, nch, s);
	copy_ui8matrix_ui8matrix(output, nrl, nrh, ncl, nch, img->Omega);
	// Save the image (debug)
	SavePGM_ui8matrix(img->Omega, nrl, nrh, ncl, nch,"output_erosion.pgm");
	free_ui8matrix(output, nrl, nrh, ncl, nch);
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