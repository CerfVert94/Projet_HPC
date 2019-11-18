/* --------------- */
/* --- morpho.h --- */
/* --------------- */
#ifndef __MORPHO_H__
#define __MORPHO_H__

#pragma message("  include  morpho.h")

typedef struct struct_elem_dim {
	//ori = origin
	long orix;
	long oriy;
	long nrow;
	long ncol;

	// th = top horizontal / bv = bottom horizontal
	// lv = left vertical / rh = right vertical
	// bord = border
	long nthbord;
	long nbhbord;
	long nlvbord;
	long nrvbord;
} struct_elem_dim;

typedef struct struct_elem_dim struct_elem_dim, *p_struct_elem_dim;
p_struct_elem_dim compute_struct_elem_dim(long orix, long oriy, long nrow, long ncol);
void free_structuring_element(p_struct_elem_dim s);

uint8** ui8matrix_dilation(uint8** input, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s);
uint8** ui8matrix_erosion(uint8** input, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s);
void image_dilation(p_image img, p_struct_elem_dim s);
void image_erosion(p_image img, p_struct_elem_dim s);

#endif /* __MORPHO_H__ */