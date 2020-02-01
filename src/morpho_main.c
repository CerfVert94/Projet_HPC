/* -------------- */
/* --- main.c --- */
/* -------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "nrdef.h"
#include "nrutil.h"
#include "img.h"
#include "mouvement.h"
#include <stdbool.h>
#include <morpho.h>
#include <mouvement.h>
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
                                        {.func_name = "ui8matrix_dilation_naive"                   , .morpho_func = ui8matrix_dilation_naive          , .pack_type = NO_PACK},
                                        {.func_name = "ui8matrix_dilation_LU3x3_O1xO1"                , .morpho_func = ui8matrix_dilation_LU3x3_O1xO1                   , .pack_type = NO_PACK},
                                        {.func_name = "ui8matrix_dilation_LU3x3_ExLU_O3"                , .morpho_func = ui8matrix_dilation_LU3x3_ExLU_O3                   , .pack_type = NO_PACK},
                                        {.func_name = "ui8matrix_dilation_LU3x3_InLU_O3"                , .morpho_func = ui8matrix_dilation_LU3x3_InLU_O3                   , .pack_type = NO_PACK},
                                        {.func_name = "ui8matrix_dilation_LU3x3_ComLU_O3"                , .morpho_func = ui8matrix_dilation_LU3x3_ComLU_O3                   , .pack_type = NO_PACK},
                                        {.func_name = "ui8matrix_dilation_LU3x3_ExLU_O3_ValAddrRR"      , .morpho_func = ui8matrix_dilation_LU3x3_ExLU_O3_ValAddrRR         , .pack_type = NO_PACK},
                                        {.func_name = "ui8matrix_dilation_LU3x3_InLU_O3_ValAddrRR"      , .morpho_func = ui8matrix_dilation_LU3x3_InLU_O3_ValAddrRR         , .pack_type = NO_PACK},
                                        {.func_name = "ui8matrix_dilation_LU3x3_ComLU_O3_ValAddrRR"      , .morpho_func = ui8matrix_dilation_LU3x3_ComLU_O3_ValAddrRR         , .pack_type = NO_PACK},
                                        {.func_name = "ui8matrix_dilation_LU3x3_InLU_O3_AddrRR"         , .morpho_func = ui8matrix_dilation_LU3x3_InLU_O3_AddrRR            , .pack_type = NO_PACK},
                                        {.func_name = "ui8matrix_dilation_LU3x3_ExLU_O3_AddrRR"         , .morpho_func = ui8matrix_dilation_LU3x3_ExLU_O3_AddrRR            , .pack_type = NO_PACK},
                                        {.func_name = "ui8matrix_dilation_LU3x3_ComLU_O3_AddrRR"         , .morpho_func = ui8matrix_dilation_LU3x3_ComLU_O3_AddrRR            , .pack_type = NO_PACK},
                                        {.func_name = "ui8matrix_dilation_LU3x3_ExLU_O3_NS"             , .morpho_func = ui8matrix_dilation_LU3x3_ExLU_O3_NS                , .pack_type = NO_PACK},
                                        {.func_name = "ui8matrix_dilation_LU3x3_ExLU_O3_RR_NS"          , .morpho_func = ui8matrix_dilation_LU3x3_ExLU_O3_RR_NS             , .pack_type = NO_PACK},
                                        {.func_name = "ui8matrix_dilation_LU3x3_InLU_O3_NS"             , .morpho_func = ui8matrix_dilation_LU3x3_InLU_O3_NS                , .pack_type = NO_PACK},
                                        {.func_name = "ui8matrix_dilation_LU3x3_InLU_O3_RR_NS"          , .morpho_func = ui8matrix_dilation_LU3x3_InLU_O3_RR_NS             , .pack_type = NO_PACK},
                                        {.func_name = "ui8matrix_dilation_LU3x3_ComLU_O3_NS"             , .morpho_func = ui8matrix_dilation_LU3x3_ComLU_O3_NS                , .pack_type = NO_PACK},
                                        {.func_name = "ui8matrix_dilation_row_pipeline"               , .morpho_func = ui8matrix_dilation_row_pipeline                  , .pack_type = NO_PACK},
                                        {.func_name = "ui8matrix_dilation_col_pipeline"               , .morpho_func = ui8matrix_dilation_col_pipeline                  , .pack_type = NO_PACK},
                                        {.func_name = "ui8matrix_dilation_col_pipeline_RR"            , .morpho_func = ui8matrix_dilation_col_pipeline_RR               , .pack_type = NO_PACK},
                                        {.func_name = "ui8matrix_dilation_pipeline_LU3x3_ExLU_O3"       , .morpho_func = ui8matrix_dilation_pipeline_LU3x3_ExLU_O3          , .pack_type = NO_PACK},
                                        {.func_name = "ui8matrix_dilation_pipeline_LU3x3_ExLU_O3_RR"    , .morpho_func = ui8matrix_dilation_pipeline_LU3x3_ExLU_O3_RR       , .pack_type = NO_PACK},
                                        {.func_name = "ui8matrix_dilation_pipeline_LU3x3_InLU_O3"       , .morpho_func = ui8matrix_dilation_pipeline_LU3x3_InLU_O3          , .pack_type = NO_PACK},
                                        {.func_name = "ui8matrix_dilation_pipeline_LU3x3_InLU_O3_RR"    , .morpho_func = ui8matrix_dilation_pipeline_LU3x3_InLU_O3_RR       , .pack_type = NO_PACK},
                                        {.func_name = "ui8matrix_dilation_pipeline2_LU3x3_ExLU_O3"             , .morpho_func = ui8matrix_dilation_pipeline2_LU3x3_ExLU_O3                , .pack_type = NO_PACK},
                                        {.func_name = "ui8matrix_dilation_pipeline2_LU3x3_ExLU_O3_RR"          , .morpho_func = ui8matrix_dilation_pipeline2_LU3x3_ExLU_O3_RR             , .pack_type = NO_PACK},
                                        {.func_name = "ui8matrix_dilation_pipeline2_LU3x3_InLU_O3"             , .morpho_func = ui8matrix_dilation_pipeline2_LU3x3_InLU_O3                , .pack_type = NO_PACK},
                                        {.func_name = "ui8matrix_dilation_pipeline2_LU3x3_InLU_O3_RR"          , .morpho_func = ui8matrix_dilation_pipeline2_LU3x3_InLU_O3_RR             , .pack_type = NO_PACK},
                                        {.func_name = "ui8matrix_dilation_divide_row_and_conquer"            , .morpho_func = ui8matrix_dilation_divide_row_and_conquer               , .pack_type = NO_PACK},
                                        {.func_name = "ui8matrix_dilation_divide_row_and_conquer_InLU_O3"    , .morpho_func = ui8matrix_dilation_divide_row_and_conquer_InLU_O3               , .pack_type = NO_PACK},
                                        {.func_name = "ui8matrix_dilation_divide_row_and_conquer_ExLU_O3"    , .morpho_func = ui8matrix_dilation_divide_row_and_conquer_ExLU_O3               , .pack_type = NO_PACK},
                                        {.func_name = "ui8matrix_dilation_divide_col_and_conquer"            , .morpho_func = ui8matrix_dilation_divide_col_and_conquer               , .pack_type = NO_PACK},
                                        {.func_name = "ui8matrix_dilation_divide_col_and_conquer_InLU_O3"    , .morpho_func = ui8matrix_dilation_divide_col_and_conquer_InLU_O3               , .pack_type = NO_PACK},
                                        {.func_name = "ui8matrix_dilation_divide_col_and_conquer_ExLU_O3"    , .morpho_func = ui8matrix_dilation_divide_col_and_conquer_ExLU_O3               , .pack_type = NO_PACK},
                                        {.func_name = "ui8matrix_dilation_LU3x3_O1xO1_OMP"            , .morpho_func = ui8matrix_dilation_LU3x3_O1xO1_OMP               , .pack_type = NO_PACK},
                                        {.func_name = "ui8matrix_dilation_LU3x3_ExLU_O3_OMP"            , .morpho_func = ui8matrix_dilation_LU3x3_ExLU_O3_OMP               , .pack_type = NO_PACK},
                                        {.func_name = "ui8matrix_dilation_LU3x3_InLU_O3_OMP"            , .morpho_func = ui8matrix_dilation_LU3x3_InLU_O3_OMP               , .pack_type = NO_PACK},
                                        {.func_name = "ui8matrix_dilation_LU3x3_ComLU_O3_OMP"            , .morpho_func = ui8matrix_dilation_LU3x3_ComLU_O3_OMP               , .pack_type = NO_PACK},
                                        {.func_name = "ui8matrix_dilation_LU3x3_InLU_O3_AddrRR_OMP"     , .morpho_func = ui8matrix_dilation_LU3x3_InLU_O3_AddrRR_OMP        , .pack_type = NO_PACK},
                                        {.func_name = "ui8matrix_dilation_LU3x3_ExLU_O3_AddrRR_OMP"     , .morpho_func = ui8matrix_dilation_LU3x3_ExLU_O3_AddrRR_OMP        , .pack_type = NO_PACK},
                                        {.func_name = "ui8matrix_dilation_LU3x3_ComLU_O3_AddrRR_OMP"     , .morpho_func = ui8matrix_dilation_LU3x3_ComLU_O3_AddrRR_OMP        , .pack_type = NO_PACK},
                                        {.func_name = "ui8matrix_dilation_pipeline_LU3x3_ExLU_O3_RR_OMP"         , .morpho_func = ui8matrix_dilation_pipeline_LU3x3_ExLU_O3_RR_OMP            , .pack_type = NO_PACK},
                                        {.func_name = "ui8matrix_dilation_pipeline2_LU3x3_InLU_O3_RR_OMP"      , .morpho_func = ui8matrix_dilation_pipeline2_LU3x3_InLU_O3_RR_OMP         , .pack_type = NO_PACK},
                                        {.func_name = "ui8matrix_dilation_divide_row_and_conquer_OMP"            , .morpho_func = ui8matrix_dilation_divide_row_and_conquer_OMP               , .pack_type = NO_PACK},
                                        {.func_name = "ui8matrix_dilation_divide_row_and_conquer_InLU_O3_OMP"    , .morpho_func = ui8matrix_dilation_divide_row_and_conquer_InLU_O3_OMP               , .pack_type = NO_PACK},
                                        {.func_name = "ui8matrix_dilation_divide_row_and_conquer_ExLU_O3_OMP"    , .morpho_func = ui8matrix_dilation_divide_row_and_conquer_ExLU_O3_OMP               , .pack_type = NO_PACK},
                                        // {.func_name = "ui8matrix_dilation_hpacked_divide_row_and_conquer"            , .morpho_func = ui8matrix_dilation_hpacked_divide_row_and_conquer               , .pack_type = HPACK},
                                        };


    
    struct morpho_set erosions[] = {
                                        {.func_name = "ui8matrix_erosion_naive"                   , .morpho_func = ui8matrix_erosion_naive          }, 
                                        {.func_name = "ui8matrix_erosion_LU3x3_O1xO1"                , .morpho_func = ui8matrix_erosion_LU3x3_O1xO1                   },
                                        {.func_name = "ui8matrix_erosion_LU3x3_ExLU_O3"                , .morpho_func = ui8matrix_erosion_LU3x3_ExLU_O3                   },
                                        {.func_name = "ui8matrix_erosion_LU3x3_InLU_O3"                , .morpho_func = ui8matrix_erosion_LU3x3_InLU_O3                   },
                                        {.func_name = "ui8matrix_erosion_LU3x3_ComLU_O3"                , .morpho_func = ui8matrix_erosion_LU3x3_ComLU_O3                   },
                                        {.func_name = "ui8matrix_erosion_LU3x3_ExLU_O3_ValAddrRR"      , .morpho_func = ui8matrix_erosion_LU3x3_ExLU_O3_ValAddrRR         },
                                        {.func_name = "ui8matrix_erosion_LU3x3_InLU_O3_ValAddrRR"      , .morpho_func = ui8matrix_erosion_LU3x3_InLU_O3_ValAddrRR         },
                                        {.func_name = "ui8matrix_erosion_LU3x3_ComLU_O3_ValAddrRR"      , .morpho_func = ui8matrix_erosion_LU3x3_ComLU_O3_ValAddrRR         },
                                        {.func_name = "ui8matrix_erosion_LU3x3_InLU_O3_AddrRR"         , .morpho_func = ui8matrix_erosion_LU3x3_InLU_O3_AddrRR            },
                                        {.func_name = "ui8matrix_erosion_LU3x3_ExLU_O3_AddrRR"         , .morpho_func = ui8matrix_erosion_LU3x3_ExLU_O3_AddrRR            },
                                        {.func_name = "ui8matrix_erosion_LU3x3_ComLU_O3_AddrRR"         , .morpho_func = ui8matrix_erosion_LU3x3_ComLU_O3_AddrRR            },
                                        {.func_name = "ui8matrix_erosion_LU3x3_ExLU_O3_NS"             , .morpho_func = ui8matrix_erosion_LU3x3_ExLU_O3_NS                },
                                        {.func_name = "ui8matrix_erosion_LU3x3_ExLU_O3_RR_NS"          , .morpho_func = ui8matrix_erosion_LU3x3_ExLU_O3_RR_NS             },
                                        {.func_name = "ui8matrix_erosion_LU3x3_InLU_O3_NS"             , .morpho_func = ui8matrix_erosion_LU3x3_InLU_O3_NS                },
                                        {.func_name = "ui8matrix_erosion_LU3x3_InLU_O3_RR_NS"          , .morpho_func = ui8matrix_erosion_LU3x3_InLU_O3_RR_NS             },
                                        {.func_name = "ui8matrix_erosion_LU3x3_ComLU_O3_NS"             , .morpho_func = ui8matrix_erosion_LU3x3_ComLU_O3_NS                },
                                        {.func_name = "ui8matrix_erosion_row_pipeline"               , .morpho_func = ui8matrix_erosion_row_pipeline                  },
                                        {.func_name = "ui8matrix_erosion_col_pipeline"               , .morpho_func = ui8matrix_erosion_col_pipeline                  },
                                        {.func_name = "ui8matrix_erosion_col_pipeline_RR"            , .morpho_func = ui8matrix_erosion_col_pipeline_RR               },
                                        {.func_name = "ui8matrix_erosion_pipeline_LU3x3_ExLU_O3"       , .morpho_func = ui8matrix_erosion_pipeline_LU3x3_ExLU_O3          },
                                        {.func_name = "ui8matrix_erosion_pipeline_LU3x3_ExLU_O3_RR"    , .morpho_func = ui8matrix_erosion_pipeline_LU3x3_ExLU_O3_RR       },
                                        {.func_name = "ui8matrix_erosion_pipeline_LU3x3_InLU_O3"       , .morpho_func = ui8matrix_erosion_pipeline_LU3x3_InLU_O3          },
                                        {.func_name = "ui8matrix_erosion_pipeline_LU3x3_InLU_O3_RR"    , .morpho_func = ui8matrix_erosion_pipeline_LU3x3_InLU_O3_RR       },
                                        {.func_name = "ui8matrix_erosion_pipeline2_LU3x3_ExLU_O3"             , .morpho_func = ui8matrix_erosion_pipeline2_LU3x3_ExLU_O3                },
                                        {.func_name = "ui8matrix_erosion_pipeline2_LU3x3_ExLU_O3_RR"          , .morpho_func = ui8matrix_erosion_pipeline2_LU3x3_ExLU_O3_RR             },
                                        {.func_name = "ui8matrix_erosion_pipeline2_LU3x3_InLU_O3"             , .morpho_func = ui8matrix_erosion_pipeline2_LU3x3_InLU_O3                },
                                        {.func_name = "ui8matrix_erosion_pipeline2_LU3x3_InLU_O3_RR"          , .morpho_func = ui8matrix_erosion_pipeline2_LU3x3_InLU_O3_RR             },
                                        {.func_name = "ui8matrix_erosion_divide_row_and_conquer"            , .morpho_func = ui8matrix_erosion_divide_row_and_conquer               },
                                        {.func_name = "ui8matrix_erosion_divide_row_and_conquer_InLU_O3"    , .morpho_func = ui8matrix_erosion_divide_row_and_conquer_InLU_O3               },
                                        {.func_name = "ui8matrix_erosion_divide_row_and_conquer_ExLU_O3"    , .morpho_func = ui8matrix_erosion_divide_row_and_conquer_ExLU_O3               },
                                        {.func_name = "ui8matrix_erosion_divide_col_and_conquer"            , .morpho_func = ui8matrix_erosion_divide_col_and_conquer               },
                                        {.func_name = "ui8matrix_erosion_divide_col_and_conquer_InLU_O3"    , .morpho_func = ui8matrix_erosion_divide_col_and_conquer_InLU_O3               },
                                        {.func_name = "ui8matrix_erosion_divide_col_and_conquer_ExLU_O3"    , .morpho_func = ui8matrix_erosion_divide_col_and_conquer_ExLU_O3               },
                                        {.func_name = "ui8matrix_erosion_LU3x3_O1xO1_OMP"            , .morpho_func = ui8matrix_erosion_LU3x3_O1xO1_OMP               },
                                        {.func_name = "ui8matrix_erosion_LU3x3_ExLU_O3_OMP"            , .morpho_func = ui8matrix_erosion_LU3x3_ExLU_O3_OMP               },
                                        {.func_name = "ui8matrix_erosion_LU3x3_InLU_O3_OMP"            , .morpho_func = ui8matrix_erosion_LU3x3_InLU_O3_OMP               },
                                        {.func_name = "ui8matrix_erosion_LU3x3_ComLU_O3_OMP"            , .morpho_func = ui8matrix_erosion_LU3x3_ComLU_O3_OMP               },
                                        {.func_name = "ui8matrix_erosion_LU3x3_InLU_O3_AddrRR_OMP"     , .morpho_func = ui8matrix_erosion_LU3x3_InLU_O3_AddrRR_OMP        },
                                        {.func_name = "ui8matrix_erosion_LU3x3_ExLU_O3_AddrRR_OMP"     , .morpho_func = ui8matrix_erosion_LU3x3_ExLU_O3_AddrRR_OMP        },
                                        {.func_name = "ui8matrix_erosion_LU3x3_ComLU_O3_AddrRR_OMP"     , .morpho_func = ui8matrix_erosion_LU3x3_ComLU_O3_AddrRR_OMP        },
                                        {.func_name = "ui8matrix_erosion_pipeline_LU3x3_ExLU_O3_RR_OMP"         , .morpho_func = ui8matrix_erosion_pipeline_LU3x3_ExLU_O3_RR_OMP            },
                                        {.func_name = "ui8matrix_erosion_pipeline2_LU3x3_InLU_O3_RR_OMP"      , .morpho_func = ui8matrix_erosion_pipeline2_LU3x3_InLU_O3_RR_OMP         },
                                        {.func_name = "ui8matrix_erosion_divide_row_and_conquer_OMP"            , .morpho_func = ui8matrix_erosion_divide_row_and_conquer_OMP               },
                                        {.func_name = "ui8matrix_erosion_divide_row_and_conquer_InLU_O3_OMP"    , .morpho_func = ui8matrix_erosion_divide_row_and_conquer_InLU_O3_OMP               },
                                        {.func_name = "ui8matrix_erosion_divide_row_and_conquer_ExLU_O3_OMP"    , .morpho_func = ui8matrix_erosion_divide_row_and_conquer_ExLU_O3_OMP               },
                                       };
    struct morpho_set sequence_sets[] = {
                                        // {.func_name = "ui8matrix_sequence_naive"                     , .morpho_func = ui8matrix_sequence_naive                        }, 
                                        // {.func_name = "ui8matrix_sequence_divide_row_and_conquer_OMP", .morpho_func = ui8matrix_sequence_divide_row_and_conquer_OMP   },
                                       };
    struct sd_set SDs_step0[] = {
                                {.func_name = "SigmaDelta_step0", .sd_func = SigmaDelta_step0_naive, .n_coeff = N, .v_min = Vmin, .v_max=Vmax},
                                };  
    struct sd_set SDs_step1[] = {
                                {.func_name = "SigmaDelta_step1", .sd_func = SigmaDelta_step1_naive, .n_coeff = N, .v_min = Vmin, .v_max=Vmax},
                                };
    struct sd_set SDs_step2[] = {
                                {.func_name = "SigmaDelta_step2", .sd_func = SigmaDelta_step2_naive, .n_coeff = N, .v_min = Vmin, .v_max =Vmax},
                                };
    struct sd_set SDs_step3[] = {
                                {.func_name = "SigmaDelta_step3", .sd_func = SigmaDelta_step3_naive, .n_coeff = N, .v_min = Vmin, .v_max =Vmax},
                                };
    struct sd_set SDs_step4[] = {
                                {.func_name = "SigmaDelta_step4", .sd_func = SigmaDelta_step4_naive, .n_coeff = N, .v_min = Vmin, .v_max =Vmax},
                                };
    struct sd_set SD_steps[5] = {
                                 {.func_name = "SigmaDelta_step0", .sd_func = SigmaDelta_step0_naive, .n_coeff = N, .v_min = Vmin, .v_max=Vmax},
                                 {.func_name = "SigmaDelta_step1", .sd_func = SigmaDelta_step1_naive, .n_coeff = N, .v_min = Vmin, .v_max=Vmax},
                                 {.func_name = "SigmaDelta_step2", .sd_func = SigmaDelta_step2_naive, .n_coeff = N, .v_min = Vmin, .v_max =Vmax},
                                 {.func_name = "SigmaDelta_step3", .sd_func = SigmaDelta_step3_naive, .n_coeff = N, .v_min = Vmin, .v_max =Vmax},
                                 {.func_name = "SigmaDelta_step4", .sd_func = SigmaDelta_step4_naive, .n_coeff = N, .v_min = Vmin, .v_max =Vmax}
                                };
    
    struct complete_sd_set completeSDs[] = {
                                    {.func_name = "SigmaDelta_naive", .sd_step0 = SigmaDelta_step0_naive, .sd_func = SigmaDelta_naive, .n_coeff = N, .v_min = Vmin, .v_max = Vmax},
                                  };

// -----------
int main(void)
// -----------
{
    // p_struct_elem_dim s5 = compute_struct_elem_dim(2,2,5,5);
    // p_struct_elem_dim s3 = compute_struct_elem_dim(1,1,3,3);
    long nb_sets;
    nb_sets = 45; 

    // test_SigmaDelta_step0("../car3/car_3000.pgm", "../car3/car_3001.pgm", SDs_step0, 1, false);
    // test_SigmaDelta_step1("../car3/car_3000.pgm", "../car3/car_3001.pgm", SDs_step1, 1, false);
    // test_SigmaDelta_step2("../car3/car_3000.pgm", "../car3/car_3001.pgm", SDs_step2, 1, false);
    // test_SigmaDelta_step3("../car3/car_3000.pgm", "../car3/car_3001.pgm", SDs_step3, 1, false);
    // test_SigmaDelta_step4("../car3/car_3000.pgm", "../car3/car_3001.pgm", SDs_step4, 1, false);
    // test_SigmaDelta("../car3/car_3000.pgm", "../car3/car_3001.pgm", completeSDs, 1, false);
    // test_erosions ("../car3/car_3000.pgm", erosions , nb_sets, false);
    // test_dilations("../car3/car_3000.pgm", dilations, nb_sets, false);
    launch_morpho_benchmark( "output/benchmark_dilation.dat", dilations  , nb_sets, 1, 1, 200, 1000, 1);
    launch_morpho_benchmark( "output/benchmark_erosion.dat" , erosions   , nb_sets, 1, 1, 200, 1000, 1);
    launch_SD_step_benchmark("output/benchmark_SD_step.dat" , SD_steps   ,       5, 1, 1, 200, 5000, 50);
    launch_SD_benchmark(     "output/benchmark_SD.dat"      , completeSDs,       1, 1, 1, 200, 5000, 50);
    
    // long ncol = 20, nrow = 20;
    // long packed_nrl, packed_nrh, packed_ncl, packed_nch, bord;
    // uint8 **a, **b, **c, **temp, **out;
    // a = ui8matrix(0, nrow, 0, ncol);
    // c = ui8matrix(0, nrow, 0, ncol);

    // for (long i = 0; i < nrow + 1; i++) {
    //     for (long j = 0; j < ncol + 1; j++) {
    //         // if (j <= i)
    //             // a[i][j] = 1;
    //         // else 
    //             a[i][j] = 0;
    //     }
    //     a[i][0] = 1;
    //     a[i][1] = 0;
    //     a[i][2] = 1;
    //     a[i][3] = 0;
    //     a[i][4] = 0;
    //     a[i][5] = 0;
    //     a[i][5] = 0;
    //     a[i][6] = 0;
    //     a[i][7] = 1;

    // }
    // b    = fcpacked_ui8matrix(0, nrow, 0, ncol, &packed_nrl, &packed_nrh, &packed_ncl, &packed_nch, &bord);
    // temp = fcpacked_ui8matrix(0, nrow, 0, ncol, &packed_nrl, &packed_nrh, &packed_ncl, &packed_nch, &bord);
    // out  = fcpacked_ui8matrix(0, nrow, 0, ncol, &packed_nrl, &packed_nrh, &packed_ncl, &packed_nch, &bord);

    // printf("%ld %ld %ld %ld\n", packed_nrl, packed_nrh, packed_ncl, packed_nch);
    // fcpack_ui8matrix_ui8matrix(a, 0, nrow, 0, ncol, packed_nrl, packed_nrh, packed_ncl, packed_nch, bord, b);
    // ui8matrix_dilation_hpacked_divide_row_and_conquer(b, packed_nrl, packed_nrh, packed_ncl, packed_nch, temp, out);
    // unfcpack_ui8matrix_ui8matrix(out, 0, nrow, 0, ncol, packed_nrl, packed_nrh, packed_ncl, packed_nch, bord, c);

    // display_ui8matrix(b, packed_nrl - bord, packed_nrh + bord, packed_ncl - bord, packed_nch + bord, "%02x", "packed");
    // display_ui8matrix(c, 0, nrow, 0, ncol, "%u", "unpacked");
    // free_packed_ui8matrix(b, packed_nrl, packed_nrh, packed_ncl, packed_nch, bord);
    // free_ui8matrix(a, 0, nrow, 0, ncol);
    // free_ui8matrix(c, 0, nrow, 0, ncol);

    // packing_test("../car3/car_3000.pgm", hpack_binary_ui8matrix, unhpack_binary_ui8matrix, "hpack", true);
    // packing_test("../car3/car_3000.pgm", vpack_binary_ui8matrix, unvpack_binary_ui8matrix, "vpack", false);

    
    // test_dilations(dilations, nb_sets, false);
    
    // launch_benchmark_compression("output/benchmark_dilation.dat", dilations, nb_sets, 1, 1, 200, 25000);
    // for (int i = -40; i < 49; i++)
    //     pack_binary_ui8matrix(NULL, i, i, i, i, &nrl, &nrh, &ncl, &nch);

    //test_erosions(erosions, nb_sets, true);
    // launch_benchmark("output/benchmark_erosion.dat", erosions, nb_sets, 200, 4000);

    // test_sequences(sequence_sets, nb_sets, false);
    // launch_benchmark("output/benchmark_sequence.dat", sequence_sets, nb_sets, 200, 2000);

    // make_testsets();
    // free_structuring_element(s3);
    // free_structuring_element(s5);

    return 0;    
}
























