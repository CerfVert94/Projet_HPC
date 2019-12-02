#ifndef  __TEST_MORPHO_H__
#define __TEST_MORPHO_H__

#pragma message("  include  test_morpho.h")

// uint8 Erosion_Test_5x5_Rect (uint8 **ppRect, p_struct_elem_dim s);
// uint8 Dilation_Test_5x5_Rect(uint8 **ppRect, p_struct_elem_dim s);
unsigned long long get_cpu_cycles(morpho_func_t morpho, uint8 **ppInput, p_struct_elem_dim s, long nrl, long nrh, long ncl, long nch, uint8 **ppOutput);
uint8 Morpho_Test_5x5_Rect(morpho_func_t morpho, uint8 **ppInput, p_struct_elem_dim s);
void test_morpho(morpho_func_t erosion, morpho_func_t dilation);
#endif /*  __TEST_MORPHO_H__ */