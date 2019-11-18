#ifndef  __TEST_MORPHO_H__
#define __TEST_MORPHO_H__

#pragma message("  include  test_morpho.h")
#define TEST_RECT_TL_CORNER     {{1,0,0,0,0},\
                                 {0,0,0,0,0},\
                                 {0,0,0,0,0},\
                                 {0,0,0,0,0},\
                                 {0,0,0,0,0}}

#define TEST_RECT_TR_CORNER     {{0,0,0,0,1},\
                                 {0,0,0,0,0},\
                                 {0,0,0,0,0},\
                                 {0,0,0,0,0},\
                                 {0,0,0,0,0}}

#define TEST_RECT_BL_CORNER     {{0,0,0,0,0},\
                                 {0,0,0,0,0},\
                                 {0,0,0,0,0},\
                                 {0,0,0,0,0},\
                                 {1,0,0,0,0}}

#define TEST_RECT_BR_CORNER     {{0,0,0,0,0},\
                                 {0,0,0,0,0},\
                                 {0,0,0,0,0},\
                                 {0,0,0,0,0},\
                                 {0,0,0,0,1}}

#define TEST_RECT_LV_EDGE       {{1,0,0,0,0},\
                                 {1,0,0,0,0},\
                                 {1,0,0,0,0},\
                                 {1,0,0,0,0},\
                                 {1,0,0,0,0}}

#define TEST_RECT_RV_EDGE       {{0,0,0,0,1},\
                                 {0,0,0,0,1},\
                                 {0,0,0,0,1},\
                                 {0,0,0,0,1},\
                                 {0,0,0,0,1}}

#define TEST_RECT_TH_EDGE       {{1,1,1,1,1},\
                                 {0,0,0,0,0},\
                                 {0,0,0,0,0},\
                                 {0,0,0,0,0},\
                                 {0,0,0,0,0}}

#define TEST_RECT_BH_EDGE       {{0,0,0,0,0},\
                                 {0,0,0,0,0},\
                                 {0,0,0,0,0},\
                                 {0,0,0,0,0},\
                                 {1,1,1,1,1}}

#define TEST_RECT_FULL          {{1,1,1,1,1},\
                                 {1,1,1,1,1},\
                                 {1,1,1,1,1},\
                                 {1,1,1,1,1},\
                                 {1,1,1,1,1}}

#define TEST_RECT_EMPTY         {{0,0,0,0,0},\
                                 {0,0,0,0,0},\
                                 {0,0,0,0,0},\
                                 {0,0,0,0,0},\
                                 {0,0,0,0,0}}

#define TEST_RECT_NONEMPTY      {{0,0,0,0,0},\
                                 {0,1,1,1,0},\
                                 {0,1,1,1,0},\
                                 {0,1,1,1,0},\
                                 {0,0,0,0,0}}
uint8 Erosion_Test_5x5_TopLeft_Corner(p_struct_elem_dim s);
uint8 Erosion_Test_5x5_TopRight_Corner(p_struct_elem_dim s);
uint8 Erosion_Test_5x5_BottomLeft_Corner(p_struct_elem_dim s);
uint8 Erosion_Test_5x5_BottomRight_Corner(p_struct_elem_dim s);
uint8 Erosion_Test_5x5_Left_Edge(p_struct_elem_dim s);
uint8 Erosion_Test_5x5_Right_Edge(p_struct_elem_dim s);
uint8 Erosion_Test_5x5_Top_Edge(p_struct_elem_dim s);
uint8 Erosion_Test_5x5_Bottom_Edge(p_struct_elem_dim s);
uint8 Erosion_Test_5x5_Empty_Rect(p_struct_elem_dim s);
uint8 Erosion_Test_5x5_NonEmpty_Rect(p_struct_elem_dim s);
uint8 Erosion_Test_5x5_Full_Rect(p_struct_elem_dim s);
uint8 Dilation_Test_5x5_TopLeft_Corner(p_struct_elem_dim s);
uint8 Dilation_Test_5x5_TopRight_Corner(p_struct_elem_dim s);
uint8 Dilation_Test_5x5_BottomLeft_Corner(p_struct_elem_dim s);
uint8 Dilation_Test_5x5_BottomRight_Corner(p_struct_elem_dim s);	
uint8 Dilation_Test_5x5_Left_Edge(p_struct_elem_dim s);
uint8 Dilation_Test_5x5_Right_Edge(p_struct_elem_dim s);
uint8 Dilation_Test_5x5_Top_Edge(p_struct_elem_dim s) ;
uint8 Dilation_Test_5x5_Bottom_Edge(p_struct_elem_dim s) ;
uint8 Dilation_Test_5x5_Empty_Rect(p_struct_elem_dim s) ;
uint8 Dilation_Test_5x5_NonEmpty_Rect(p_struct_elem_dim s) ;
uint8 Dilation_Test_5x5_Full_Rect(p_struct_elem_dim s);
#endif /*  __TEST_MORPHO_H__ */