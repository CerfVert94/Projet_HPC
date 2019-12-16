/* --------------- */
/* --- morpho.h --- */
/* --------------- */
#ifndef __MORPHO_H__
#define __MORPHO_H__

#pragma message("  include  morpho.h")

typedef struct struct_elem_dim {
	//ori = origin
	long x0;
	long y0;
	long nrow;
	long ncol;

	long nrl;
	long nrh;
	long ncl;
	long nch;
} struct_elem_dim;

typedef struct struct_elem_dim struct_elem_dim, *p_struct_elem_dim;
p_struct_elem_dim compute_struct_elem_dim(long x0, long y0, long nrow, long ncol);
void free_structuring_element(p_struct_elem_dim s);


typedef void (*morpho_func_t)(uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput);

void ui8matrix_dilation_naive(uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput);
void ui8matrix_dilation_LU3x3_O1xO1 (uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput);
void ui8matrix_dilation_LU3x3_O3xO1 (uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput);
void ui8matrix_dilation_LU3x3_O3xO1_RR (uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput);
void ui8matrix_dilation_LU3x3_O3xO3 (uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput);
void ui8matrix_dilation_LU3x3_O3xO3_RR (uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput);
void ui8matrix_dilation_pipeline3x3 (uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput);
void ui8matrix_dilation_pipeline_LU3x3_O3xO1 (uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput);
void ui8matrix_dilation_pipeline_LU3x3_O3xO1_RR (uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput);

void ui8matrix_dilation_LU5x5 (uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput);


void ui8matrix_erosion_naive (uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput);
void ui8matrix_erosion_LU3x3 (uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput);
void ui8matrix_erosion_LU5x5 (uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput);
void ui8matrix_erosion_LU3x3_O3 (uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput);
void ui8matrix_erosion_LU3x3_O3_RR (uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput);
void ui8matrix_erosion_LU5x5_O3_RR (uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput);


void ui8matrix_sequence_naive(uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput);
void ui8matrix_sequence_naive_inline(uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput);
void ui8matrix_sequence_LU3x3(uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput);





/******************/
/*** DEPRECATED ***/
/******************/

// uint8  erosion_naive(uint8** ppInput, long row, long col, p_struct_elem_dim s);
// uint8 dilation_naive(uint8** ppInput, long row, long col, p_struct_elem_dim s);



// void ui8matrix_lambda_morpho(morpho_func_t morpho, long order, uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput);
// void dilation_naive(uint8** ppInput, long row, long col, p_struct_elem_dim s, uint8 **ppOutput);
// void  erosion_naive(uint8** ppInput, long row, long col, p_struct_elem_dim s, uint8 **ppOutput);


void image_chain_processing(p_image img, p_struct_elem_dim s, int idx);
#endif /* __MORPHO_H__ */