#include "nrdef.h"
#include "nrutil.h"
#include <malloc.h>

// Without optimisation
p_struct_elem create_structuring_element(long orix, long oriy, long nrow, long ncol)
{
	p_struct_elem s;

	s = (p_struct_elem) malloc(sizeof(struct_elem));

	s->m = ui8matrix(0, nrow, 0, ncol);

	s->nrow = nrow;
	s->ncol = ncol;
	s->orix= orix;
	s->oriy= oriy;
}

uint8 **erosion(p_image img, p_struct_elem s) 
{
    
}

uint8 **dilation(uint8** input, uint8** mask, long orix, long oriy) 
{
    
}