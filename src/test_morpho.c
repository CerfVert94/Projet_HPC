#include <nrdef.h>
#include <nrutil.h>
#include <img.h>
#include <morpho.h>
#include <stdio.h>
#include <stdlib.h>
#include <test_morpho.h>
#include <util.h>

const char * nom_func;
void test_morpho()
{
	const long TEST_SE_NROW     = 5L;
	const long TEST_SE_NCOL     = 5L;
	const long TEST_SE_ORIGIN_X = 2L;
	const long TEST_SE_ORIGIN_Y = 2L;
	p_struct_elem_dim s = compute_struct_elem_dim(TEST_SE_ORIGIN_X, TEST_SE_ORIGIN_Y, 
	 											  TEST_SE_NROW    , TEST_SE_NCOL);

	
	// image_erosion(create_image("morpho_test.pgm"), s);
	// image_dilation(create_image("morpho_test.pgm"), s);
	puts("Checking the structuring element used for unit tests.");

	if (s->nrow != 5 || s->ncol != 5) {
		fprintf(stderr, "\tIncorrect dimension : \n\t\trequires (%ld, %ld) but has (%ld, %ld)\n", TEST_SE_NROW, TEST_SE_NCOL, s->nrow, s->ncol);
		exit(EXIT_FAILURE);
	}

	if (s->oriy != 2 || s->orix != 2) {
		fprintf(stderr, "\tIncorrect origin  : \n\t\trequires (%ld, %ld) but has (%ld, %ld)\n", TEST_SE_ORIGIN_X, TEST_SE_ORIGIN_Y, s->nrow, s->ncol);
		exit(EXIT_FAILURE);
	}
	puts("The structuring element for unit tests has the correct dimension.");
	
	UNIT_TEST( Erosion_Test_5x5_TopLeft_Corner(s)	   == 0   , nom_func, "");
	UNIT_TEST( Erosion_Test_5x5_TopRight_Corner(s)	   == 0   , nom_func, "");
	UNIT_TEST( Erosion_Test_5x5_BottomLeft_Corner(s)   == 0   , nom_func, "");
	UNIT_TEST( Erosion_Test_5x5_BottomRight_Corner(s)  == 0   , nom_func, "");
	UNIT_TEST( Erosion_Test_5x5_Left_Edge(s)           == 0   , nom_func, "");
	UNIT_TEST( Erosion_Test_5x5_Right_Edge(s)          == 0   , nom_func, "");
	UNIT_TEST( Erosion_Test_5x5_Top_Edge(s)            == 0   , nom_func, "");
	UNIT_TEST( Erosion_Test_5x5_Bottom_Edge(s)         == 0   , nom_func, "");
	UNIT_TEST( Erosion_Test_5x5_Empty_Rect(s)          == 0   , nom_func, "");
	UNIT_TEST( Erosion_Test_5x5_NonEmpty_Rect(s)       == 0   , nom_func, "");
	UNIT_TEST( Erosion_Test_5x5_Full_Rect(s)           == 1   , nom_func, "");

	UNIT_TEST(Dilation_Test_5x5_TopLeft_Corner(s)      == 1   , nom_func, "");
	UNIT_TEST(Dilation_Test_5x5_TopRight_Corner(s)     == 1   , nom_func, "");
	UNIT_TEST(Dilation_Test_5x5_BottomLeft_Corner(s)   == 1   , nom_func, "");
	UNIT_TEST(Dilation_Test_5x5_BottomRight_Corner(s)  == 1   , nom_func, "");
	UNIT_TEST(Dilation_Test_5x5_Left_Edge(s)           == 1   , nom_func, "");
	UNIT_TEST(Dilation_Test_5x5_Right_Edge(s)          == 1   , nom_func, "");
	UNIT_TEST(Dilation_Test_5x5_Top_Edge(s)            == 1   , nom_func, "");
	UNIT_TEST(Dilation_Test_5x5_Bottom_Edge(s)         == 1   , nom_func, "");
	UNIT_TEST(Dilation_Test_5x5_NonEmpty_Rect(s)       == 1   , nom_func, "");
	UNIT_TEST(Dilation_Test_5x5_Full_Rect(s)           == 1   , nom_func, "");
	UNIT_TEST(Dilation_Test_5x5_Empty_Rect(s)          == 0   , nom_func, "");
	puts("Morpho : Passed all tests.");
}

#define INIT_TEST_RECT(TEST_RECT)					uint8 rect[5][5] = TEST_RECT;\
													uint8 *pRows[5] = {rect[0], rect[1], rect[2], rect[3], rect[4]};\
													uint8 **input = pRows, **output;\
													uint8 pixel

#define MEMORIZE_FUNC_CALL 							nom_func = __func__

#define TEST_MORPHO(func, input, output, pixel, s)	output = func(input, 0, 4, 0, 4, s);\
													pixel = output[s->oriy][s->orix];\
													free_ui8matrix(output, 0, 4, 0, 4)\

uint8 Erosion_Test_5x5_TopLeft_Corner(p_struct_elem_dim s)  {
	INIT_TEST_RECT(TEST_RECT_TL_CORNER);
	MEMORIZE_FUNC_CALL;
	TEST_MORPHO(ui8matrix_erosion, input, output, pixel, s);
	return pixel;
}
uint8 Erosion_Test_5x5_TopRight_Corner(p_struct_elem_dim s) {
	INIT_TEST_RECT(TEST_RECT_TR_CORNER);
	MEMORIZE_FUNC_CALL;
	TEST_MORPHO(ui8matrix_erosion, input, output, pixel, s);
	return pixel;
}
uint8 Erosion_Test_5x5_BottomLeft_Corner(p_struct_elem_dim s) {
	INIT_TEST_RECT(TEST_RECT_BL_CORNER);
	MEMORIZE_FUNC_CALL;
	TEST_MORPHO(ui8matrix_erosion, input, output, pixel, s);
	return pixel;
}
uint8 Erosion_Test_5x5_BottomRight_Corner(p_struct_elem_dim s) {
	INIT_TEST_RECT(TEST_RECT_BR_CORNER);
	MEMORIZE_FUNC_CALL;
	TEST_MORPHO(ui8matrix_erosion, input, output, pixel, s);
	return pixel;
}
uint8 Erosion_Test_5x5_Left_Edge(p_struct_elem_dim s) {
	INIT_TEST_RECT(TEST_RECT_LV_EDGE);
	MEMORIZE_FUNC_CALL;
	TEST_MORPHO(ui8matrix_erosion, input, output, pixel, s);
	return pixel;
}
uint8 Erosion_Test_5x5_Right_Edge(p_struct_elem_dim s) {
	INIT_TEST_RECT(TEST_RECT_RV_EDGE);
	MEMORIZE_FUNC_CALL;
	TEST_MORPHO(ui8matrix_erosion, input, output, pixel, s);
	return pixel;
}
uint8 Erosion_Test_5x5_Top_Edge(p_struct_elem_dim s) {
	INIT_TEST_RECT(TEST_RECT_TH_EDGE);
	MEMORIZE_FUNC_CALL;
	TEST_MORPHO(ui8matrix_erosion, input, output, pixel, s);
	return pixel;
}
uint8 Erosion_Test_5x5_Bottom_Edge(p_struct_elem_dim s) {
	INIT_TEST_RECT(TEST_RECT_TH_EDGE);
	MEMORIZE_FUNC_CALL;
	TEST_MORPHO(ui8matrix_erosion, input, output, pixel, s);
	return pixel;
}
uint8 Erosion_Test_5x5_Empty_Rect(p_struct_elem_dim s) {
	INIT_TEST_RECT(TEST_RECT_EMPTY);
	MEMORIZE_FUNC_CALL;
	TEST_MORPHO(ui8matrix_erosion, input, output, pixel, s);
	return pixel;
}
uint8 Erosion_Test_5x5_NonEmpty_Rect(p_struct_elem_dim s) {
	INIT_TEST_RECT(TEST_RECT_NONEMPTY);
	MEMORIZE_FUNC_CALL;
	TEST_MORPHO(ui8matrix_erosion, input, output, pixel, s);
	return pixel;
}
uint8 Erosion_Test_5x5_Full_Rect(p_struct_elem_dim s) {
	INIT_TEST_RECT(TEST_RECT_FULL);
	MEMORIZE_FUNC_CALL;
	TEST_MORPHO(ui8matrix_erosion, input, output, pixel, s);
	return pixel;
}
uint8 Dilation_Test_5x5_TopLeft_Corner(p_struct_elem_dim s) 
{
	INIT_TEST_RECT(TEST_RECT_TL_CORNER);
	MEMORIZE_FUNC_CALL;
	TEST_MORPHO(ui8matrix_dilation, input, output, pixel, s);
	return pixel;
}
uint8 Dilation_Test_5x5_TopRight_Corner(p_struct_elem_dim s) {
	INIT_TEST_RECT(TEST_RECT_TR_CORNER);
	MEMORIZE_FUNC_CALL;
	TEST_MORPHO(ui8matrix_dilation, input, output, pixel, s);
	return pixel;
}
uint8 Dilation_Test_5x5_BottomLeft_Corner(p_struct_elem_dim s) {
	INIT_TEST_RECT(TEST_RECT_BL_CORNER);
	MEMORIZE_FUNC_CALL;
	TEST_MORPHO(ui8matrix_dilation, input, output, pixel, s);
	return pixel;
}
uint8 Dilation_Test_5x5_BottomRight_Corner(p_struct_elem_dim s) {
	INIT_TEST_RECT(TEST_RECT_BR_CORNER);
	MEMORIZE_FUNC_CALL;
	TEST_MORPHO(ui8matrix_dilation, input, output, pixel, s);
	return pixel;
}
uint8 Dilation_Test_5x5_Left_Edge(p_struct_elem_dim s) {
	INIT_TEST_RECT(TEST_RECT_LV_EDGE);
	MEMORIZE_FUNC_CALL;
	TEST_MORPHO(ui8matrix_dilation, input, output, pixel, s);
	return pixel;
}
uint8 Dilation_Test_5x5_Right_Edge(p_struct_elem_dim s) {
	INIT_TEST_RECT(TEST_RECT_RV_EDGE);
	MEMORIZE_FUNC_CALL;
	TEST_MORPHO(ui8matrix_dilation, input, output, pixel, s);
	return pixel;
}
uint8 Dilation_Test_5x5_Top_Edge(p_struct_elem_dim s) {
	INIT_TEST_RECT(TEST_RECT_TH_EDGE);
	MEMORIZE_FUNC_CALL;
	TEST_MORPHO(ui8matrix_dilation, input, output, pixel, s);
	return pixel;
}
uint8 Dilation_Test_5x5_Bottom_Edge(p_struct_elem_dim s) {
	INIT_TEST_RECT(TEST_RECT_BH_EDGE);
	MEMORIZE_FUNC_CALL;
	TEST_MORPHO(ui8matrix_dilation, input, output, pixel, s);
	return pixel;
}
uint8 Dilation_Test_5x5_Empty_Rect(p_struct_elem_dim s) {
	INIT_TEST_RECT(TEST_RECT_EMPTY);
	MEMORIZE_FUNC_CALL;
	TEST_MORPHO(ui8matrix_dilation, input, output, pixel, s);
	return pixel;
}
uint8 Dilation_Test_5x5_NonEmpty_Rect(p_struct_elem_dim s) {
	INIT_TEST_RECT(TEST_RECT_NONEMPTY);
	MEMORIZE_FUNC_CALL;
	TEST_MORPHO(ui8matrix_dilation, input, output, pixel, s);
	return pixel;
}
uint8 Dilation_Test_5x5_Full_Rect(p_struct_elem_dim s) {
	INIT_TEST_RECT(TEST_RECT_FULL);
	MEMORIZE_FUNC_CALL;
	TEST_MORPHO(ui8matrix_dilation, input, output, pixel, s);
	return pixel;
}

