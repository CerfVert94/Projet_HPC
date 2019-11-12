#include "nrdef.h"
#include "nrutil.h"
#include <img.h>
#include <malloc.h>
#include <morpho.h>


	uint8 mask3x3_plus[3][3] = {{0,1,0},
								{1,1,1},
								{0,1,0}};


	uint8 mask3x3_ltri[3][3] = {{1,0,0},
								{1,1,0},
								{1,1,1}};

	uint8 mask3x3_htri[3][3] = {{1,1,1},
								{0,1,1},
								{0,0,0}};

	uint8 mask3x3_flat[3][3] = {{1,1,1},
								{1,1,1},
								{1,1,1}};

	uint8 mask2x2_ltri[2][2] = {{1,0},
								{1,1}};							

	uint8 mask2x2_htri[2][2] = {{1,1},
								{0,1}};							

	uint8 mask2x2_flat[2][2] = {{1,1},
								{0,1}};							
		

p_struct_elem create_structuring_element(long orix, long oriy, long nrow, long ncol)
{
	p_struct_elem s;

	s = (p_struct_elem) malloc(sizeof(struct_elem));

	// Create a nrow x ncol uint8 matrix.
	s->m = ui8matrix(0, nrow - 1, 0, ncol - 1);

	// Define dimension.
	s->nrow = nrow;
	s->ncol = ncol;
	// Define origin.
	s->orix= orix;
	s->oriy= oriy;
}

// Without optimisation
uint8 **erosion(p_image img, p_struct_elem s) 
{
	// bord = border
	// tv = top vertical / bv = bottom vertical
	// lh = left horizontal / rh = right horizontal
    long ntvbord, nlhbord, nbvbord, nrhbord;
	uint8 **input, **output;

	// Compute the size of borders.
	ntvbord = s->oriy;
	nbvbord = (s->nrow - 1) - s->oriy;
	nlhbord = s->orix;
	nrhbord = (s->ncol - 1) - s->orix;
	
	// Create input and output matrix.
	// Add edges to the input matrix.
	input = ui8matrix(img->nrl - ntvbord, img->nrh + nbvbord, img->ncl - nlhbord, img->nch + nrhbord);
	printf("%ld => %ld /  %ld => %ld\n",img->nrl - ntvbord, img->nrh + nbvbord, img->ncl - nlhbord, img->nch + nrhbord);
	output = ui8matrix(img->nrl, img->nrh, img->ncl, img->nch);
	// input[-1][-1] = 1;
	// input[8][8] = 1;
	// input[-1][8] = 1;
	// input[8][-1] = 1;
	display_ui8matrix(input, img->nrl - ntvbord, img->nrh + nbvbord, img->ncl - nlhbord, img->nch + nrhbord, "%03u ", "erosion_input_matrix");

	// Copy the original image to the input matrix.
	copy_ui8matrix_ui8matrix(img->I, img->nrl, img->nrh, img->ncl, img->nch, input);	
	display_ui8matrix(input, img->nrl - ntvbord, img->nrh + nbvbord, img->ncl - nlhbord, img->nch + nrhbord, "%03u ", "copy");


}

uint8 **dilation(uint8** input, uint8** mask, long orix, long oriy) 
{
    
}
void test_morpho()
{		
	p_struct_elem s = create_structuring_element(2,0,3,3);
	for (int i = 0; i < s->nrow; i++) {
		for (int j = 0; j < s->ncol; j++) {
			s->m[i][j] = mask3x3_plus[i][j];
		}
	}
	display_ui8matrix(s->m, 0, s->nrow - 1, 0, s->ncol - 1, "%u ", "mask3x3_plus");
	

	erosion(create_image("../car3/8x8.pgm"), s);
}