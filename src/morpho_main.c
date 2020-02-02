/* -------------- */
/* --- main.c --- */
/* -------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>

#include "nrdef.h"
#include "vnrdef.h"
#include "nrutil.h"
#include "mynrutil.h"
#include "vnrutil.h"
#include "myvnrutil.h"
#include "util.h"
#include "img.h"
#include "img_SIMD.h"

#include <morpho.h>
#include <mouvement.h>
#include <mouvement_SIMD.h>
//#include "mymacro.h"
//#include "test_simd1.h"
#include <test_mouvement.h>
#include <test_morpho.h>
#include <util.h>
#include <benchmark.h>

// ============
void info(void)
// ============
{
#ifdef ENABLE_BENCHMARK
    puts("#############################");
    puts("mode Benchmark ON & DEBUG OFF");
    puts("#############################");
#else
    puts("#############################");
    puts("mode Benchmark OFF & DEBUG ON");
    puts("#############################");
#endif
}

void movement_detection()
{
    // char filename[128];
    // // p_struct_elem_dim s = compute_struct_elem_dim(1,1,3,3);
    // int idx = 3000;
    // for (int i = 0; i < 200 - 1; i++) 
    // {
	//     sprintf(filename, "../car3/car_%d.pgm", idx + i);
    //     printf("../car3/car_%d.pgm\n", idx + i);
    //     p_image t_1 = create_image(filename);
	//     sprintf(filename, "../car3/car_%d.pgm", idx + i + 1);
	//     p_image t = create_image(filename);
    //     // SigmaDelta_step0(t);
    //     // SigmaDelta_step0(t_1);
    //     SigmaDelta(t, t_1);
	//     sprintf(filename, "../md/car_%d.pgm", idx + i + 1);
    //     // ui8matrix_dilation_divide_row_and_conquer_OMP(t->E, t->nrl + BORD, t->nrh - BORD, t->ncl + BORD, t->nch - BORD, s, t->Omega);
    //     // ui8matrix_sequence_LU3x3(t->E, t->nrl + BORD, t->nrh - BORD, t->ncl + BORD, t->nch - BORD, s, t->Omega);
    //     binary_to_octal_ui8matrix(t->Omega, t->nrl + BORD, t->nrh - BORD, t->ncl + BORD, t->nch - BORD);
    //     SavePGM_ui8matrix(t->Omega, t->nrl + BORD, t->nrh - BORD, t->ncl + BORD, t->nch - BORD, filename);
    //     free_image(t);
    //     free_image(t_1);        
    // }
    //info();
    // free_structuring_element(s);
}


void make_testsets(){
    // p_struct_elem_dim s5 = compute_struct_elem_dim(2,2,5,5);
    // p_struct_elem_dim s3 = compute_struct_elem_dim(1,1,3,3);
    long nrl, nrh, ncl, nch, bord, x, y;
    long size  = 500, noise_size = size;
    bord = 2;
    uint8** image = LoadPGM_ui8matrix("bwimg_r1.pgm", &nrl, &nrh, &ncl, &nch);
    uint8** ppInput  = ui8matrix(nrl - bord, nrh + bord, ncl - bord, nch + bord);
    uint8** ppOutput = ui8matrix(nrl - bord, nrh + bord, ncl - bord, nch + bord);
    memset_ui8matrix(ppInput, 0, nrl - bord, nrh + bord, ncl - bord, nch + bord);
    memset_ui8matrix(ppOutput, 0, nrl - bord, nrh + bord, ncl - bord, nch + bord);

    // copy_ui8matrix_ui8matrix(image, nrl, nrh, ncl, nch, ppInput);
    // octal_to_binary_ui8matrix(ppInput, nrl, nrh, ncl, nch);
    // display_ui8matrix(ppInput, nrl - bord, nrh + bord, ncl - bord, nch + bord, "%u,", "uint8 morpho5x5_test_input[] = {\\");
    // printf("};");

    // ui8matrix_erosion_naive(ppInput, nrl, nrh, ncl, nch, ppOutput);
    // display_ui8matrix(ppOutput,  nrl - bord, nrh + bord, ncl - bord, nch + bord, "%u,", "uint8 erosion5x5_test_output[] = {\\");
    // printf("};");
    // // getchar();
    // ui8matrix_dilation_naive(ppInput, nrl, nrh, ncl, nch, ppOutput);
    // display_ui8matrix(ppOutput,  nrl - bord, nrh + bord, ncl - bord, nch + bord, "%u,", "uint8 dilation5x5_test_output[] = {\\");
    // printf("};");
    // free_ui8matrix(ppInput, nrl - bord, nrh + bord, ncl - bord, nch + bord);
    // free_ui8matrix(ppOutput, nrl - bord, nrh + bord, ncl - bord, nch + bord);
    
    bord = 1;
    ppInput  = ui8matrix(nrl - bord, nrh + bord, ncl - bord, nch + bord);
    ppOutput = ui8matrix(nrl - bord, nrh + bord, ncl - bord, nch + bord);
    memset_ui8matrix(ppInput, 0, nrl - bord, nrh + bord, ncl - bord, nch + bord);
    memset_ui8matrix(ppOutput, 0, nrl - bord, nrh + bord, ncl - bord, nch + bord);

    copy_ui8matrix_ui8matrix(image, nrl, nrh, ncl, nch, ppInput);
    // octal_to_binary_ui8matrix(ppInput, nrl, nrh, ncl, nch);
    display_ui8matrix(ppInput, nrl - bord, nrh + bord, ncl - bord, nch + bord, "%u,", "uint8 morpho3x3_r2_test_input[] = {\\");
    printf("};\n");

    ui8matrix_erosion_naive(ppInput, nrl, nrh, ncl, nch, NULL, ppOutput);
    display_ui8matrix(ppOutput,  nrl - bord, nrh + bord, ncl - bord, nch + bord, "%u,", "uint8 erosion3x3_r2_test_output[] = {\\");
    printf("};\n");
    // getchar();
    ui8matrix_dilation_naive(ppInput, nrl, nrh, ncl, nch, NULL,ppOutput);
    display_ui8matrix(ppOutput,  nrl - bord, nrh + bord, ncl - bord, nch + bord, "%u,", "uint8 dilation3x3_r2_test_output[] = {\\");
    printf("};\n");
    free_ui8matrix(ppInput, nrl - bord, nrh + bord, ncl - bord, nch + bord);
    free_ui8matrix(ppOutput, nrl - bord, nrh + bord, ncl - bord, nch + bord);
    // free_structuring_element(s3);
    // free_structuring_element(s5);
}

struct morpho_set dilations[] = {
                                        // {.func_name = "ui8matrix_dilation_naive"                   , .morpho_func = ui8matrix_dilation_naive                     , .pack_type = NO_PACK, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_LU3x3_O1xO1"                , .morpho_func = ui8matrix_dilation_LU3x3_O1xO1                   , .pack_type = NO_PACK, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_LU3x3_ExLU_O3"                , .morpho_func = ui8matrix_dilation_LU3x3_ExLU_O3                   , .pack_type = NO_PACK, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_LU3x3_InLU_O3"                , .morpho_func = ui8matrix_dilation_LU3x3_InLU_O3                   , .pack_type = NO_PACK, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_LU3x3_ComLU_O3"                , .morpho_func = ui8matrix_dilation_LU3x3_ComLU_O3                   , .pack_type = NO_PACK, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_LU3x3_ExLU_O3_ValAddrRR"      , .morpho_func = ui8matrix_dilation_LU3x3_ExLU_O3_ValAddrRR         , .pack_type = NO_PACK, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_LU3x3_InLU_O3_ValAddrRR"      , .morpho_func = ui8matrix_dilation_LU3x3_InLU_O3_ValAddrRR         , .pack_type = NO_PACK, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_LU3x3_ComLU_O3_ValAddrRR"      , .morpho_func = ui8matrix_dilation_LU3x3_ComLU_O3_ValAddrRR         , .pack_type = NO_PACK, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_LU3x3_InLU_O3_AddrRR"         , .morpho_func = ui8matrix_dilation_LU3x3_InLU_O3_AddrRR            , .pack_type = NO_PACK, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_LU3x3_ExLU_O3_AddrRR"         , .morpho_func = ui8matrix_dilation_LU3x3_ExLU_O3_AddrRR            , .pack_type = NO_PACK, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_LU3x3_ComLU_O3_AddrRR"         , .morpho_func = ui8matrix_dilation_LU3x3_ComLU_O3_AddrRR            , .pack_type = NO_PACK, .instr_type = SCALAR},
                                        // {.func_name = "ui8matrix_dilation_LU3x3_ExLU_O3_NS"             , .morpho_func = ui8matrix_dilation_LU3x3_ExLU_O3_NS                , .pack_type = NO_PACK, .instr_type = SCALAR},
                                        // {.func_name = "ui8matrix_dilation_LU3x3_ExLU_O3_RR_NS"          , .morpho_func = ui8matrix_dilation_LU3x3_ExLU_O3_RR_NS             , .pack_type = NO_PACK, .instr_type = SCALAR},
                                        // {.func_name = "ui8matrix_dilation_LU3x3_InLU_O3_NS"             , .morpho_func = ui8matrix_dilation_LU3x3_InLU_O3_NS                , .pack_type = NO_PACK, .instr_type = SCALAR},
                                        // {.func_name = "ui8matrix_dilation_LU3x3_InLU_O3_RR_NS"          , .morpho_func = ui8matrix_dilation_LU3x3_InLU_O3_RR_NS             , .pack_type = NO_PACK, .instr_type = SCALAR},
                                        // {.func_name = "ui8matrix_dilation_LU3x3_ComLU_O3_NS"             , .morpho_func = ui8matrix_dilation_LU3x3_ComLU_O3_NS                , .pack_type = NO_PACK, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_row_pipeline"               , .morpho_func = ui8matrix_dilation_row_pipeline                  , .pack_type = NO_PACK, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_col_pipeline"               , .morpho_func = ui8matrix_dilation_col_pipeline                  , .pack_type = NO_PACK, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_col_pipeline_RR"            , .morpho_func = ui8matrix_dilation_col_pipeline_RR               , .pack_type = NO_PACK, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_pipeline_LU3x3_ExLU_O3"       , .morpho_func = ui8matrix_dilation_pipeline_LU3x3_ExLU_O3          , .pack_type = NO_PACK, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_pipeline_LU3x3_ExLU_O3_RR"    , .morpho_func = ui8matrix_dilation_pipeline_LU3x3_ExLU_O3_RR       , .pack_type = NO_PACK, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_pipeline_LU3x3_InLU_O3"       , .morpho_func = ui8matrix_dilation_pipeline_LU3x3_InLU_O3          , .pack_type = NO_PACK, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_pipeline_LU3x3_InLU_O3_RR"    , .morpho_func = ui8matrix_dilation_pipeline_LU3x3_InLU_O3_RR       , .pack_type = NO_PACK, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_pipeline2_LU3x3_ExLU_O3"             , .morpho_func = ui8matrix_dilation_pipeline2_LU3x3_ExLU_O3                , .pack_type = NO_PACK, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_pipeline2_LU3x3_ExLU_O3_RR"          , .morpho_func = ui8matrix_dilation_pipeline2_LU3x3_ExLU_O3_RR             , .pack_type = NO_PACK, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_pipeline2_LU3x3_InLU_O3"             , .morpho_func = ui8matrix_dilation_pipeline2_LU3x3_InLU_O3                , .pack_type = NO_PACK, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_pipeline2_LU3x3_InLU_O3_RR"          , .morpho_func = ui8matrix_dilation_pipeline2_LU3x3_InLU_O3_RR             , .pack_type = NO_PACK, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_divide_row_and_conquer"            , .morpho_func = ui8matrix_dilation_divide_row_and_conquer               , .pack_type = NO_PACK, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_divide_row_and_conquer_InLU_O3"    , .morpho_func = ui8matrix_dilation_divide_row_and_conquer_InLU_O3               , .pack_type = NO_PACK, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_divide_row_and_conquer_ExLU_O3"    , .morpho_func = ui8matrix_dilation_divide_row_and_conquer_ExLU_O3               , .pack_type = NO_PACK, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_divide_col_and_conquer"            , .morpho_func = ui8matrix_dilation_divide_col_and_conquer               , .pack_type = NO_PACK, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_divide_col_and_conquer_InLU_O3"    , .morpho_func = ui8matrix_dilation_divide_col_and_conquer_InLU_O3               , .pack_type = NO_PACK, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_divide_col_and_conquer_ExLU_O3"    , .morpho_func = ui8matrix_dilation_divide_col_and_conquer_ExLU_O3               , .pack_type = NO_PACK, .instr_type = SCALAR},
                                        // {.func_name = "ui8matrix_dilation_LU3x3_O1xO1_OMP"            , .morpho_func = ui8matrix_dilation_LU3x3_O1xO1_OMP               , .pack_type = NO_PACK, .instr_type = SCALAR},
                                        // {.func_name = "ui8matrix_dilation_LU3x3_ExLU_O3_OMP"            , .morpho_func = ui8matrix_dilation_LU3x3_ExLU_O3_OMP               , .pack_type = NO_PACK, .instr_type = SCALAR},
                                        // {.func_name = "ui8matrix_dilation_LU3x3_InLU_O3_OMP"            , .morpho_func = ui8matrix_dilation_LU3x3_InLU_O3_OMP               , .pack_type = NO_PACK, .instr_type = SCALAR},
                                        // {.func_name = "ui8matrix_dilation_LU3x3_ComLU_O3_OMP"            , .morpho_func = ui8matrix_dilation_LU3x3_ComLU_O3_OMP               , .pack_type = NO_PACK, .instr_type = SCALAR},
                                        // {.func_name = "ui8matrix_dilation_LU3x3_InLU_O3_AddrRR_OMP"     , .morpho_func = ui8matrix_dilation_LU3x3_InLU_O3_AddrRR_OMP        , .pack_type = NO_PACK, .instr_type = SCALAR},
                                        // {.func_name = "ui8matrix_dilation_LU3x3_ExLU_O3_AddrRR_OMP"     , .morpho_func = ui8matrix_dilation_LU3x3_ExLU_O3_AddrRR_OMP        , .pack_type = NO_PACK, .instr_type = SCALAR},
                                        // {.func_name = "ui8matrix_dilation_LU3x3_ComLU_O3_AddrRR_OMP"     , .morpho_func = ui8matrix_dilation_LU3x3_ComLU_O3_AddrRR_OMP        , .pack_type = NO_PACK, .instr_type = SCALAR},
                                        // {.func_name = "ui8matrix_dilation_pipeline_LU3x3_ExLU_O3_RR_OMP"         , .morpho_func = ui8matrix_dilation_pipeline_LU3x3_ExLU_O3_RR_OMP            , .pack_type = NO_PACK, .instr_type = SCALAR},
                                        // {.func_name = "ui8matrix_dilation_pipeline2_LU3x3_InLU_O3_RR_OMP"      , .morpho_func = ui8matrix_dilation_pipeline2_LU3x3_InLU_O3_RR_OMP         , .pack_type = NO_PACK, .instr_type = SCALAR},
                                        // {.func_name = "ui8matrix_dilation_divide_row_and_conquer_OMP"            , .morpho_func = ui8matrix_dilation_divide_row_and_conquer_OMP               , .pack_type = NO_PACK, .instr_type = SCALAR},
                                        // {.func_name = "ui8matrix_dilation_divide_row_and_conquer_InLU_O3_OMP"    , .morpho_func = ui8matrix_dilation_divide_row_and_conquer_InLU_O3_OMP               , .pack_type = NO_PACK, .instr_type = SCALAR},
                                        // {.func_name = "ui8matrix_dilation_divide_row_and_conquer_ExLU_O3_OMP"    , .morpho_func = ui8matrix_dilation_divide_row_and_conquer_ExLU_O3_OMP               , .pack_type = NO_PACK, .instr_type = SCALAR},
                                        // {.func_name = "ui8matrix_dilation_hpacked_divide_row_and_conquer"            , .morpho_func = ui8matrix_dilation_hpacked_divide_row_and_conquer               , .pack_type = HPACK, .instr_type = SCALAR},
                                        };


    
    struct morpho_set erosions[] = {
                                        // {.func_name = "ui8matrix_erosion_naive"                   , .morpho_func = ui8matrix_erosion_naive          },.instr_type = SCALAR, 
                                        {.func_name = "ui8matrix_erosion_LU3x3_O1xO1"                , .morpho_func = ui8matrix_erosion_LU3x3_O1xO1                   ,.instr_type = SCALAR},
                                        {.func_name = "ui8matrix_erosion_LU3x3_ExLU_O3"                , .morpho_func = ui8matrix_erosion_LU3x3_ExLU_O3                   ,.instr_type = SCALAR},
                                        {.func_name = "ui8matrix_erosion_LU3x3_InLU_O3"                , .morpho_func = ui8matrix_erosion_LU3x3_InLU_O3                   ,.instr_type = SCALAR},
                                        {.func_name = "ui8matrix_erosion_LU3x3_ComLU_O3"                , .morpho_func = ui8matrix_erosion_LU3x3_ComLU_O3                   ,.instr_type = SCALAR},
                                        {.func_name = "ui8matrix_erosion_LU3x3_ExLU_O3_ValAddrRR"      , .morpho_func = ui8matrix_erosion_LU3x3_ExLU_O3_ValAddrRR         ,.instr_type = SCALAR},
                                        {.func_name = "ui8matrix_erosion_LU3x3_InLU_O3_ValAddrRR"      , .morpho_func = ui8matrix_erosion_LU3x3_InLU_O3_ValAddrRR         ,.instr_type = SCALAR},
                                        {.func_name = "ui8matrix_erosion_LU3x3_ComLU_O3_ValAddrRR"      , .morpho_func = ui8matrix_erosion_LU3x3_ComLU_O3_ValAddrRR         ,.instr_type = SCALAR},
                                        {.func_name = "ui8matrix_erosion_LU3x3_InLU_O3_AddrRR"         , .morpho_func = ui8matrix_erosion_LU3x3_InLU_O3_AddrRR            ,.instr_type = SCALAR},
                                        {.func_name = "ui8matrix_erosion_LU3x3_ExLU_O3_AddrRR"         , .morpho_func = ui8matrix_erosion_LU3x3_ExLU_O3_AddrRR            ,.instr_type = SCALAR},
                                        {.func_name = "ui8matrix_erosion_LU3x3_ComLU_O3_AddrRR"         , .morpho_func = ui8matrix_erosion_LU3x3_ComLU_O3_AddrRR            ,.instr_type = SCALAR},
                                        {.func_name = "ui8matrix_erosion_LU3x3_ExLU_O3_NS"             , .morpho_func = ui8matrix_erosion_LU3x3_ExLU_O3_NS                ,.instr_type = SCALAR},
                                        {.func_name = "ui8matrix_erosion_LU3x3_ExLU_O3_RR_NS"          , .morpho_func = ui8matrix_erosion_LU3x3_ExLU_O3_RR_NS             ,.instr_type = SCALAR},
                                        {.func_name = "ui8matrix_erosion_LU3x3_InLU_O3_NS"             , .morpho_func = ui8matrix_erosion_LU3x3_InLU_O3_NS                ,.instr_type = SCALAR},
                                        {.func_name = "ui8matrix_erosion_LU3x3_InLU_O3_RR_NS"          , .morpho_func = ui8matrix_erosion_LU3x3_InLU_O3_RR_NS             ,.instr_type = SCALAR},
                                        {.func_name = "ui8matrix_erosion_LU3x3_ComLU_O3_NS"             , .morpho_func = ui8matrix_erosion_LU3x3_ComLU_O3_NS                ,.instr_type = SCALAR},
                                        {.func_name = "ui8matrix_erosion_row_pipeline"               , .morpho_func = ui8matrix_erosion_row_pipeline                  ,.instr_type = SCALAR},
                                        {.func_name = "ui8matrix_erosion_col_pipeline"               , .morpho_func = ui8matrix_erosion_col_pipeline                  ,.instr_type = SCALAR},
                                        {.func_name = "ui8matrix_erosion_col_pipeline_RR"            , .morpho_func = ui8matrix_erosion_col_pipeline_RR               ,.instr_type = SCALAR},
                                        {.func_name = "ui8matrix_erosion_pipeline_LU3x3_ExLU_O3"       , .morpho_func = ui8matrix_erosion_pipeline_LU3x3_ExLU_O3          ,.instr_type = SCALAR},
                                        {.func_name = "ui8matrix_erosion_pipeline_LU3x3_ExLU_O3_RR"    , .morpho_func = ui8matrix_erosion_pipeline_LU3x3_ExLU_O3_RR       ,.instr_type = SCALAR},
                                        {.func_name = "ui8matrix_erosion_pipeline_LU3x3_InLU_O3"       , .morpho_func = ui8matrix_erosion_pipeline_LU3x3_InLU_O3          ,.instr_type = SCALAR},
                                        {.func_name = "ui8matrix_erosion_pipeline_LU3x3_InLU_O3_RR"    , .morpho_func = ui8matrix_erosion_pipeline_LU3x3_InLU_O3_RR       ,.instr_type = SCALAR},
                                        {.func_name = "ui8matrix_erosion_pipeline2_LU3x3_ExLU_O3"             , .morpho_func = ui8matrix_erosion_pipeline2_LU3x3_ExLU_O3                ,.instr_type = SCALAR},
                                        {.func_name = "ui8matrix_erosion_pipeline2_LU3x3_ExLU_O3_RR"          , .morpho_func = ui8matrix_erosion_pipeline2_LU3x3_ExLU_O3_RR             ,.instr_type = SCALAR},
                                        {.func_name = "ui8matrix_erosion_pipeline2_LU3x3_InLU_O3"             , .morpho_func = ui8matrix_erosion_pipeline2_LU3x3_InLU_O3                ,.instr_type = SCALAR},
                                        {.func_name = "ui8matrix_erosion_pipeline2_LU3x3_InLU_O3_RR"          , .morpho_func = ui8matrix_erosion_pipeline2_LU3x3_InLU_O3_RR             ,.instr_type = SCALAR},
                                        {.func_name = "ui8matrix_erosion_divide_row_and_conquer"            , .morpho_func = ui8matrix_erosion_divide_row_and_conquer               ,.instr_type = SCALAR},
                                        {.func_name = "ui8matrix_erosion_divide_row_and_conquer_InLU_O3"    , .morpho_func = ui8matrix_erosion_divide_row_and_conquer_InLU_O3               ,.instr_type = SCALAR},
                                        {.func_name = "ui8matrix_erosion_divide_row_and_conquer_ExLU_O3"    , .morpho_func = ui8matrix_erosion_divide_row_and_conquer_ExLU_O3               ,.instr_type = SCALAR},
                                        {.func_name = "ui8matrix_erosion_divide_col_and_conquer"            , .morpho_func = ui8matrix_erosion_divide_col_and_conquer               ,.instr_type = SCALAR},
                                        {.func_name = "ui8matrix_erosion_divide_col_and_conquer_InLU_O3"    , .morpho_func = ui8matrix_erosion_divide_col_and_conquer_InLU_O3               ,.instr_type = SCALAR},
                                        {.func_name = "ui8matrix_erosion_divide_col_and_conquer_ExLU_O3"    , .morpho_func = ui8matrix_erosion_divide_col_and_conquer_ExLU_O3               ,.instr_type = SCALAR},
                                        {.func_name = "ui8matrix_erosion_LU3x3_O1xO1_OMP"            , .morpho_func = ui8matrix_erosion_LU3x3_O1xO1_OMP               ,.instr_type = SCALAR},
                                        {.func_name = "ui8matrix_erosion_LU3x3_ExLU_O3_OMP"            , .morpho_func = ui8matrix_erosion_LU3x3_ExLU_O3_OMP               ,.instr_type = SCALAR},
                                        {.func_name = "ui8matrix_erosion_LU3x3_InLU_O3_OMP"            , .morpho_func = ui8matrix_erosion_LU3x3_InLU_O3_OMP               ,.instr_type = SCALAR},
                                        {.func_name = "ui8matrix_erosion_LU3x3_ComLU_O3_OMP"            , .morpho_func = ui8matrix_erosion_LU3x3_ComLU_O3_OMP               ,.instr_type = SCALAR},
                                        {.func_name = "ui8matrix_erosion_LU3x3_InLU_O3_AddrRR_OMP"     , .morpho_func = ui8matrix_erosion_LU3x3_InLU_O3_AddrRR_OMP        ,.instr_type = SCALAR},
                                        {.func_name = "ui8matrix_erosion_LU3x3_ExLU_O3_AddrRR_OMP"     , .morpho_func = ui8matrix_erosion_LU3x3_ExLU_O3_AddrRR_OMP        ,.instr_type = SCALAR},
                                        {.func_name = "ui8matrix_erosion_LU3x3_ComLU_O3_AddrRR_OMP"     , .morpho_func = ui8matrix_erosion_LU3x3_ComLU_O3_AddrRR_OMP        ,.instr_type = SCALAR},
                                        {.func_name = "ui8matrix_erosion_pipeline_LU3x3_ExLU_O3_RR_OMP"         , .morpho_func = ui8matrix_erosion_pipeline_LU3x3_ExLU_O3_RR_OMP            ,.instr_type = SCALAR},
                                        {.func_name = "ui8matrix_erosion_pipeline2_LU3x3_InLU_O3_RR_OMP"      , .morpho_func = ui8matrix_erosion_pipeline2_LU3x3_InLU_O3_RR_OMP         ,.instr_type = SCALAR},
                                        {.func_name = "ui8matrix_erosion_divide_row_and_conquer_OMP"            , .morpho_func = ui8matrix_erosion_divide_row_and_conquer_OMP               ,.instr_type = SCALAR},
                                        {.func_name = "ui8matrix_erosion_divide_row_and_conquer_InLU_O3_OMP"    , .morpho_func = ui8matrix_erosion_divide_row_and_conquer_InLU_O3_OMP               ,.instr_type = SCALAR},
                                        {.func_name = "ui8matrix_erosion_divide_row_and_conquer_ExLU_O3_OMP"    , .morpho_func = ui8matrix_erosion_divide_row_and_conquer_ExLU_O3_OMP               ,.instr_type = SCALAR},
                                       };
    struct morpho_set sequence_sets[] = {
                                        // {.func_name = "ui8matrix_sequence_naive"                     , .morpho_func = ui8matrix_sequence_naive                        , .instr_type = SCALAR}, 
                                        // {.func_name = "ui8matrix_sequence_divide_row_and_conquer_OMP", .morpho_func = ui8matrix_sequence_divide_row_and_conquer_OMP   , .instr_type = SCALAR},
                                       };
    struct sd_set SDs_step0[] = {
                                {.func_name = "SigmaDelta_step0_naive", .sd_func = SigmaDelta_step0_naive, .n_coeff = N, .v_min = Vmin, .v_max=Vmax, .instr_type = SCALAR},
                                {.func_name = "SigmaDelta_step0_mem", .sd_func = SigmaDelta_step0_mem, .n_coeff = N, .v_min = Vmin, .v_max=Vmax, .instr_type = SCALAR},
                                {.func_name = "SigmaDelta_step0_SIMD", .vec_sd_func = SigmaDelta_step0_SIMD, .n_coeff = N, .v_min = Vmin, .v_max=Vmax, .instr_type = SIMD},
                                }; 






    struct sd_set SDs_step1[] = {
                                {.func_name = "SigmaDelta_step1_naive"  , .sd_func = SigmaDelta_step1_naive, .n_coeff = N, .v_min = Vmin, .v_max=Vmax, .instr_type = SCALAR},
                                {/*Best :0.22*/.func_name = "SigmaDelta_step1_InLU_O3_NoIf", .sd_func = SigmaDelta_step1_InLU_O3_NoIf, .n_coeff = N, .v_min = Vmin, .v_max=Vmax, .instr_type = SCALAR},
                                {              .func_name = "SigmaDelta_step1_ExLU_O3_NoIf", .sd_func = SigmaDelta_step1_ExLU_O3_NoIf, .n_coeff = N, .v_min = Vmin, .v_max=Vmax, .instr_type = SCALAR},
                                {.func_name = "SigmaDelta_step1_SIMD", .vec_sd_func = SigmaDelta_step1_SIMD, .n_coeff = N, .v_min = Vmin, .v_max=Vmax, .instr_type = SIMD},
                                };
    struct sd_set SDs_step2[] = {
                                {.func_name = "SigmaDelta_step2_naive", .sd_func = SigmaDelta_step2_naive, .n_coeff = N, .v_min = Vmin, .v_max =Vmax, .instr_type = SCALAR},
                                {/*Best :0.33*/.func_name = "SigmaDelta_step2_InLU_O3_bitop", .sd_func = SigmaDelta_step2_InLU_O3_bitop, .n_coeff = N, .v_min = Vmin, .v_max=Vmax, .instr_type = SCALAR},
                                {              .func_name = "SigmaDelta_step2_ExLU_O3_bitop", .sd_func = SigmaDelta_step2_ExLU_O3_bitop, .n_coeff = N, .v_min = Vmin, .v_max=Vmax, .instr_type = SCALAR},
                                {.func_name = "SigmaDelta_step2_SIMD", .vec_sd_func = SigmaDelta_step2_SIMD, .n_coeff = N, .v_min = Vmin, .v_max=Vmax, .instr_type = SIMD},
                                }
                                ;
    struct sd_set SDs_step3[] = {
                                {.func_name = "SigmaDelta_step3_naive", .sd_func = SigmaDelta_step3_naive, .n_coeff = N, .v_min = Vmin, .v_max =Vmax, .instr_type = SCALAR},
                                {/*Best :2.94*/.func_name = "SigmaDelta_step3_InLU_O3_NoIf", .sd_func = SigmaDelta_step3_InLU_O3_NoIf, .n_coeff = N, .v_min = Vmin, .v_max =Vmax, .instr_type = SCALAR},
                                {              .func_name = "SigmaDelta_step3_ExLU_O3_NoIf", .sd_func = SigmaDelta_step3_ExLU_O3_NoIf, .n_coeff = N, .v_min = Vmin, .v_max =Vmax, .instr_type = SCALAR},
                                {.func_name = "SigmaDelta_step3_SIMD", .vec_sd_func = SigmaDelta_step3_SIMD, .n_coeff = N, .v_min = Vmin, .v_max=Vmax, .instr_type = SIMD},
                                };
    struct sd_set SDs_step4[] = {
                                {.func_name = "SigmaDelta_step4_naive", .sd_func = SigmaDelta_step4_naive, .n_coeff = N, .v_min = Vmin, .v_max =Vmax, .instr_type = SCALAR},
                                {/*Best :0.18*/.func_name = "SigmaDelta_step4_InLU_O3_NoIf", .sd_func = SigmaDelta_step4_InLU_O3_NoIf, .n_coeff = N, .v_min = Vmin, .v_max =Vmax, .instr_type = SCALAR},
                                {              .func_name = "SigmaDelta_step4_ExLU_O3_NoIf", .sd_func = SigmaDelta_step4_ExLU_O3_NoIf, .n_coeff = N, .v_min = Vmin, .v_max =Vmax, .instr_type = SCALAR},
                                {.func_name = "SigmaDelta_step4_SIMD", .vec_sd_func = SigmaDelta_step4_SIMD, .n_coeff = N, .v_min = Vmin, .v_max=Vmax, .instr_type = SIMD},
                                };
                                
                                
    struct sd_set SD_steps[] =  {
                                 {.func_name = "SigmaDelta_step0_naive", .sd_func = SigmaDelta_step0_naive, .n_coeff = N, .v_min = Vmin, .v_max=Vmax ,.instr_type = SCALAR},
                                 {.func_name = "SigmaDelta_step1_naive", .sd_func = SigmaDelta_step1_naive, .n_coeff = N, .v_min = Vmin, .v_max=Vmax ,.instr_type = SCALAR},
                                 {.func_name = "SigmaDelta_step2_naive", .sd_func = SigmaDelta_step2_naive, .n_coeff = N, .v_min = Vmin, .v_max =Vmax ,.instr_type = SCALAR},
                                 {.func_name = "SigmaDelta_step3_naive", .sd_func = SigmaDelta_step3_naive, .n_coeff = N, .v_min = Vmin, .v_max =Vmax ,.instr_type = SCALAR},
                                 {.func_name = "SigmaDelta_step4_naive", .sd_func = SigmaDelta_step4_naive, .n_coeff = N, .v_min = Vmin, .v_max =Vmax ,.instr_type = SCALAR},
                                 {.func_name = "SigmaDelta_step0_mem", .sd_func = SigmaDelta_step0_mem, .n_coeff = N, .v_min = Vmin, .v_max=Vmax ,.instr_type = SCALAR},
                                 {/*Best :0.22*/.func_name = "SigmaDelta_step1_InLU_O3_NoIf", .sd_func = SigmaDelta_step1_InLU_O3_NoIf, .n_coeff = N, .v_min = Vmin, .v_max=Vmax ,.instr_type = SCALAR},
                                 {              .func_name = "SigmaDelta_step1_ExLU_O3_NoIf", .sd_func = SigmaDelta_step1_ExLU_O3_NoIf, .n_coeff = N, .v_min = Vmin, .v_max=Vmax ,.instr_type = SCALAR},
                                 {/*Best :0.33*/.func_name = "SigmaDelta_step2_InLU_O3_bitop", .sd_func = SigmaDelta_step2_InLU_O3_bitop, .n_coeff = N, .v_min = Vmin, .v_max=Vmax ,.instr_type = SCALAR},
                                 {              .func_name = "SigmaDelta_step2_ExLU_O3_bitop", .sd_func = SigmaDelta_step2_ExLU_O3_bitop, .n_coeff = N, .v_min = Vmin, .v_max=Vmax ,.instr_type = SCALAR},
                                 {/*Best :2.94*/.func_name = "SigmaDelta_step3_InLU_O3_NoIf", .sd_func = SigmaDelta_step3_InLU_O3_NoIf, .n_coeff = N, .v_min = Vmin, .v_max =Vmax ,.instr_type = SCALAR},
                                 {              .func_name = "SigmaDelta_step3_ExLU_O3_NoIf", .sd_func = SigmaDelta_step3_ExLU_O3_NoIf, .n_coeff = N, .v_min = Vmin, .v_max =Vmax ,.instr_type = SCALAR},
                                 {/*Best :0.18*/.func_name = "SigmaDelta_step4_InLU_O3_NoIf", .sd_func = SigmaDelta_step4_InLU_O3_NoIf, .n_coeff = N, .v_min = Vmin, .v_max =Vmax ,.instr_type = SCALAR},
                                 {              .func_name = "SigmaDelta_step4_ExLU_O3_NoIf", .sd_func = SigmaDelta_step4_ExLU_O3_NoIf, .n_coeff = N, .v_min = Vmin, .v_max =Vmax ,.instr_type = SCALAR},
                                 {.func_name = "SigmaDelta_step0_SIMD", .vec_sd_func = SigmaDelta_step0_SIMD, .n_coeff = N, .v_min = Vmin, .v_max=Vmax, .instr_type = SIMD},
                                 {.func_name = "SigmaDelta_step1_SIMD", .vec_sd_func = SigmaDelta_step1_SIMD, .n_coeff = N, .v_min = Vmin, .v_max=Vmax, .instr_type = SIMD},
                                 {.func_name = "SigmaDelta_step2_SIMD", .vec_sd_func = SigmaDelta_step2_SIMD, .n_coeff = N, .v_min = Vmin, .v_max=Vmax, .instr_type = SIMD},
                                 {.func_name = "SigmaDelta_step3_SIMD", .vec_sd_func = SigmaDelta_step3_SIMD, .n_coeff = N, .v_min = Vmin, .v_max=Vmax, .instr_type = SIMD},
                                 {.func_name = "SigmaDelta_step4_SIMD", .vec_sd_func = SigmaDelta_step4_SIMD, .n_coeff = N, .v_min = Vmin, .v_max=Vmax, .instr_type = SIMD},
                                };
    
    struct complete_sd_set completeSDs[] = {
                                    {.func_name = "SigmaDelta_naive", .sd_step0 = SigmaDelta_step0_naive, .sd_func = SigmaDelta_naive, .n_coeff = N, .v_min = Vmin, .v_max = Vmax, .instr_type = SCALAR},                                  
                                    {.func_name = "SigmaDelta_best", .sd_step0 = SigmaDelta_step0_mem, .sd_func = SigmaDelta_best, .n_coeff = N, .v_min = Vmin, .v_max = Vmax, .instr_type = SCALAR},
                                  };

// -----------
int main(void)
// -----------
{
    long nrl, nrh;
    int v0, v1;
    vuint8 **vI_t0 = LoadPGM_vui8matrix("../car3/car_3000.pgm", &nrl, &nrh, &v0, &v1);
    vuint8 **vM_t0 = vui8matrix(0,0,0,0);
    vuint8 **vM_t1 = vui8matrix(0,0,0,0);

    long nb_sets;
    nb_sets = 27; 
    // display_vui8matrix(vI_t0, 0, 0, 0, 0, "%4u", "vimage");
    // SigmaDelta_step0_SIMD(vM_t0, vI_t0, vM_t1, 0,0,0,0, N, Vmin, Vmax);
    // display_vui8matrix(vM_t0, 0, 0, 0, 0, "%4u", "vM_t0");
    // display_vui8matrix(vI_t0, 0, 0, 0, 0, "%4u", "vI_t0");
    // display_vui8matrix(vM_t1, 0, 0, 0, 0, "%4u", "vM_t1");
    test_SigmaDelta_step0("../car3/car_3000.pgm", "../car3/car_3001.pgm", SDs_step0, 3, false);
    // launch_SD_step_benchmark("output/benchmark_sdstep0.dat"             , SDs_step0, 3, 1, 1, 200, 5000, 10);
    // test_SigmaDelta_step1("../car3/car_3000.pgm", "../car3/car_3001.pgm", SDs_step1, 4, false);
    // launch_SD_step_benchmark("output/benchmark_sdstep1.dat"             , SDs_step1, 4, 1, 1, 200, 5000, 10);
    // test_SigmaDelta_step2("../car3/car_3000.pgm", "../car3/car_3001.pgm", SDs_step2, 4, false);
    // launch_SD_step_benchmark("output/benchmark_sdstep2.dat"             , SDs_step2, 4, 1, 1, 200, 5000, 10);
    // test_SigmaDelta_step3("../car3/car_3000.pgm", "../car3/car_3001.pgm", SDs_step3, 4, false);
    // launch_SD_step_benchmark("output/benchmark_sdstep3.dat"             , SDs_step3, 4, 1, 1, 200, 5000, 10);
    // test_SigmaDelta_step4("../car3/car_3000.pgm", "../car3/car_3001.pgm", SDs_step4, 4, false);
    // launch_SD_step_benchmark("output/benchmark_sdstep4.dat"             , SDs_step4, 4, 1, 1, 200, 5000, 10);
    // test_SigmaDelta("../car3/car_3000.pgm", "../car3/car_3001.pgm", completeSDs, 2, false);
    // launch_SD_step_benchmark("output/benchmark_SD_step.dat"       , SD_steps   ,       14, 1, 1, 200, 5000, 10);
    // launch_SD_benchmark(     "output/benchmark_SD.dat"            , completeSDs,       2, 1, 1, 200, 5000, 10);
    // test_erosions ("../car3/car_3000.pgm", erosions , nb_sets, false);
    // test_dilations("../car3/car_3000.pgm", dilations, nb_sets, false);
    // launch_morpho_benchmark( "output/benchmark_dilation.dat", dilations  , nb_sets, 1, 1, 10, 1000, 10);
    // launch_morpho_benchmark( "output/benchmark_erosion.dat" , erosions   , nb_sets, 1, 1, 200, 1000, 1);

    return 0;    
}
























