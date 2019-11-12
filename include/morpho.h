/* --------------- */
/* --- morpho.h --- */
/* --------------- */
#ifndef __MORPHO_H__
#define __MORPHO_H__

#pragma message("  include  morpho.h")

typedef struct struct_elem {
	long orix;
	long oriy;
	long nrow;
	long ncol;

	// uint8** m;
} struct_elem;


#define COPY_MASK2STRUCT_ELEM(mask, s) copy_ui8matrix_ui8matrix(mask, 0, s->nrow, 0, s->ncol, s->m)
typedef struct struct_elem struct_elem, *p_struct_elem;
void test_morpho();

#endif /* __MORPHO_H__ */