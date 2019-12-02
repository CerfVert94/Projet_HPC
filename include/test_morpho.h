#ifndef  __TEST_MORPHO_H__
#define __TEST_MORPHO_H__

#pragma message("  include  test_morpho.h")

uint8 Erosion_Test_5x5_Rect (uint8 **ppRect, p_struct_elem_dim s);
uint8 Dilation_Test_5x5_Rect(uint8 **ppRect, p_struct_elem_dim s);

void test_morpho();
#endif /*  __TEST_MORPHO_H__ */