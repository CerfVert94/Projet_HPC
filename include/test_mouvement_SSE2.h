/* -------------------- */
/* --- test_mouvement_SSE2.h --- */
/* -------------------- */

#ifndef __TEST_MOUVEMENT_SSE2_H__
#define __TEST_MOUVEMENT_SSE2_H__


#pragma message("  include  test_morpho.h")

// A struct type that contains a morpho function and a structuring element
struct SD_SSE_set{
    char func_name[128];
    void (*morpho_func)(uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput);
    struct struct_elem_dim *s;
};

#endif /* __TEST_MOUVEMENT_SSE2_H__ */
