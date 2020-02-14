/* -------------- */
/* --- main.c --- */
/* -------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>
#include <emmintrin.h>
#include <tmmintrin.h>

#include "nrdef.h"
#include "vnrdef.h"
#include "nrutil.h"
#include "mynrutil.h"
#include "vnrutil.h"
#include "myvnrutil.h"
#include "util.h"
#include "img.h"
#include "img_SIMD.h"
#include "simd_macro.h"

#include <morpho.h>
#include <morpho_SIMD.h>
#include <mouvement.h>
#include <mouvement_SIMD.h>
#include <test_mouvement.h>
#include <test_morpho.h>
#include <util.h>
#include <benchmark.h>


struct morpho_set dilations[] = {
                                        {.func_name = "ui8matrix_dilation_naive"                         , .morpho_func = ui8matrix_dilation_naive                     , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation5_naive"                        , .morpho_func = ui8matrix_dilation5_naive                     , .pack_type = NO_PACK, .op_type = FUSION, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_SIMD_naive"                    , .vec_morpho_func = ui8matrix_dilation_SIMD_naive                     , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SIMD},
                                        {.func_name = "ui8matrix_dilation_SIMD_RR_row"                   , .vec_morpho_func = ui8matrix_dilation_SIMD_RR_row                     , .pack_type = NO_PACK,  .op_type = NORMAL, .instr_type = SIMD},
                                        {.func_name = "ui8matrix_dilation_SIMD_InLU_O3_AddrRR"                   , .vec_morpho_func = ui8matrix_dilation_SIMD_InLU_O3_AddrRR                     , .pack_type = NO_PACK,  .op_type = NORMAL, .instr_type = SIMD},
                                        {.func_name = "ui8matrix_dilation5_SIMD_InLU_O3_ValAddrRR"      , .vec_morpho_func = ui8matrix_dilation5_SIMD_InLU_O3_ValAddrRR                     , .pack_type = NO_PACK,  .op_type = FUSION, .instr_type = SIMD},
                                        {.func_name = "ui8matrix_dilation_SIMD_InLU_O3_ValAddrRR"       , .vec_morpho_func = ui8matrix_dilation_SIMD_InLU_O3_ValAddrRR                     , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SIMD},
                                        {.func_name = "ui8matrix_dilation_dilation_SIMD_FO"              , .vec_morpho_func = ui8matrix_dilation_dilation_SIMD_FO                     , .pack_type = NO_PACK,  .op_type = FUSION, .instr_type = SIMD},
                                        {.func_name = "ui8matrix_dilation_dilation_SIMD_FO_RR_row"       , .vec_morpho_func = ui8matrix_dilation_dilation_SIMD_FO_RR_row                     , .pack_type = NO_PACK, .op_type = FUSION, .instr_type = SIMD}, 
                                        {.func_name = "ui8matrix_dilation_LU3x3_O1xO1"                   , .morpho_func = ui8matrix_dilation_LU3x3_O1xO1                   , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_LU3x3_ExLU_O3"                 , .morpho_func = ui8matrix_dilation_LU3x3_ExLU_O3                   , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_LU3x3_InLU_O3"                 , .morpho_func = ui8matrix_dilation_LU3x3_InLU_O3                   , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_LU3x3_ComLU_O3"                , .morpho_func = ui8matrix_dilation_LU3x3_ComLU_O3                   , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_LU3x3_ExLU_O3_ValAddrRR"       , .morpho_func = ui8matrix_dilation_LU3x3_ExLU_O3_ValAddrRR         , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_LU3x3_InLU_O3_ValAddrRR"       , .morpho_func = ui8matrix_dilation_LU3x3_InLU_O3_ValAddrRR         , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_LU3x3_ComLU_O3_ValAddrRR"      , .morpho_func = ui8matrix_dilation_LU3x3_ComLU_O3_ValAddrRR         , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_LU3x3_InLU_O3_AddrRR"          , .morpho_func = ui8matrix_dilation_LU3x3_InLU_O3_AddrRR            , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_LU3x3_ExLU_O3_AddrRR"          , .morpho_func = ui8matrix_dilation_LU3x3_ExLU_O3_AddrRR            , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_LU3x3_ComLU_O3_AddrRR"         , .morpho_func = ui8matrix_dilation_LU3x3_ComLU_O3_AddrRR            , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_LU3x3_ExLU_O3_NS"              , .morpho_func = ui8matrix_dilation_LU3x3_ExLU_O3_NS                , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_LU3x3_ExLU_O3_RR_NS"           , .morpho_func = ui8matrix_dilation_LU3x3_ExLU_O3_RR_NS             , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_LU3x3_InLU_O3_NS"              , .morpho_func = ui8matrix_dilation_LU3x3_InLU_O3_NS                , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_LU3x3_InLU_O3_RR_NS"           , .morpho_func = ui8matrix_dilation_LU3x3_InLU_O3_RR_NS             , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_LU3x3_ComLU_O3_NS"             , .morpho_func = ui8matrix_dilation_LU3x3_ComLU_O3_NS                , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_row_pipeline"                  , .morpho_func = ui8matrix_dilation_row_pipeline                  , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_col_pipeline"                  , .morpho_func = ui8matrix_dilation_col_pipeline                  , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_col_pipeline_RR"               , .morpho_func = ui8matrix_dilation_col_pipeline_RR               , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_pipeline_LU3x3_ExLU_O3"        , .morpho_func = ui8matrix_dilation_pipeline_LU3x3_ExLU_O3          , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_pipeline_LU3x3_ExLU_O3_RR"     , .morpho_func = ui8matrix_dilation_pipeline_LU3x3_ExLU_O3_RR       , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_pipeline_LU3x3_InLU_O3"        , .morpho_func = ui8matrix_dilation_pipeline_LU3x3_InLU_O3          , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_pipeline_LU3x3_InLU_O3_RR"     , .morpho_func = ui8matrix_dilation_pipeline_LU3x3_InLU_O3_RR       , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_pipeline2_LU3x3_ExLU_O3"       , .morpho_func = ui8matrix_dilation_pipeline2_LU3x3_ExLU_O3                , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_pipeline2_LU3x3_ExLU_O3_RR"    , .morpho_func = ui8matrix_dilation_pipeline2_LU3x3_ExLU_O3_RR             , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_pipeline2_LU3x3_InLU_O3"       , .morpho_func = ui8matrix_dilation_pipeline2_LU3x3_InLU_O3                , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_pipeline2_LU3x3_InLU_O3_RR"    , .morpho_func = ui8matrix_dilation_pipeline2_LU3x3_InLU_O3_RR             , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_row_so"        , .morpho_func = ui8matrix_dilation_row_so               , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_row_so_InLU_O3", .morpho_func = ui8matrix_dilation_row_so_InLU_O3               , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_row_so_ExLU_O3", .morpho_func = ui8matrix_dilation_row_so_ExLU_O3               , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_col_so"        , .morpho_func = ui8matrix_dilation_col_so               , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_col_so_InLU_O3", .morpho_func = ui8matrix_dilation_col_so_InLU_O3               , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_col_so_ExLU_O3", .morpho_func = ui8matrix_dilation_col_so_ExLU_O3               , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_LU3x3_O1xO1_OMP"            , .morpho_func = ui8matrix_dilation_LU3x3_O1xO1_OMP               , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_LU3x3_ExLU_O3_OMP"            , .morpho_func = ui8matrix_dilation_LU3x3_ExLU_O3_OMP               , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_LU3x3_InLU_O3_OMP"            , .morpho_func = ui8matrix_dilation_LU3x3_InLU_O3_OMP               , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_LU3x3_ComLU_O3_OMP"            , .morpho_func = ui8matrix_dilation_LU3x3_ComLU_O3_OMP               , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_LU3x3_InLU_O3_AddrRR_OMP"     , .morpho_func = ui8matrix_dilation_LU3x3_InLU_O3_AddrRR_OMP        , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_LU3x3_ExLU_O3_AddrRR_OMP"     , .morpho_func = ui8matrix_dilation_LU3x3_ExLU_O3_AddrRR_OMP        , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_LU3x3_ComLU_O3_AddrRR_OMP"     , .morpho_func = ui8matrix_dilation_LU3x3_ComLU_O3_AddrRR_OMP        , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_pipeline_LU3x3_ExLU_O3_RR_OMP"         , .morpho_func = ui8matrix_dilation_pipeline_LU3x3_ExLU_O3_RR_OMP            , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_pipeline2_LU3x3_InLU_O3_RR_OMP"      , .morpho_func = ui8matrix_dilation_pipeline2_LU3x3_InLU_O3_RR_OMP         , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_row_so_OMP"            , .morpho_func = ui8matrix_dilation_row_so_OMP               , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_row_so_InLU_O3_OMP"    , .morpho_func = ui8matrix_dilation_row_so_InLU_O3_OMP               , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR},
                                        {.func_name = "ui8matrix_dilation_row_so_ExLU_O3_OMP"    , .morpho_func = ui8matrix_dilation_row_so_ExLU_O3_OMP               , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR},
                                        //{.func_name = "ui8matrix_dilation_hpacked_row_so"            , .morpho_func = ui8matrix_dilation_hpacked_row_so               , .pack_type = HPACK, .instr_type = SCALAR},
                                        };


    
    struct morpho_set erosions[] = { 
                
		                        {.func_name = "ui8matrix_erosion_SIMD_row_so"  , .vec_morpho_func = ui8matrix_erosion_SIMD_row_so     , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SIMD},
                                        {.func_name = "ui8matrix_erosion_SIMD_col_so"  , .vec_morpho_func = ui8matrix_erosion_SIMD_col_so     , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SIMD},
                                        {.func_name = "ui8matrix_erosion_SIMD_pipeline2_LU3x3_InLU_O3_RR"  , .vec_morpho_func = ui8matrix_erosion_SIMD_pipeline2_LU3x3_InLU_O3_RR     , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SIMD},
                                        {.func_name = "ui8matrix_erosion_SIMD_InLU_O3_AddrRR"  , .vec_morpho_func = ui8matrix_erosion_SIMD_InLU_O3_AddrRR     , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SIMD},
                                        {.func_name = "ui8matrix_erosion_SIMD_InLU_O3_ValAddrRR"  , .vec_morpho_func = ui8matrix_erosion_SIMD_InLU_O3_ValAddrRR     , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SIMD},
                                        {.func_name = "ui8matrix_erosion_SIMD_col_pipeline"   , .vec_morpho_func = ui8matrix_erosion_SIMD_col_pipeline     , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SIMD},
                                        {.func_name = "ui8matrix_erosion_SIMD_col_pipeline_RR"   , .vec_morpho_func = ui8matrix_erosion_SIMD_col_pipeline_RR     , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SIMD},
                                        {.func_name = "ui8matrix_erosion_SIMD_naive"                   , .vec_morpho_func = ui8matrix_erosion_SIMD_naive                      , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SIMD},
                                        {.func_name = "ui8matrix_erosion_SIMD_RR_row"                   , .vec_morpho_func = ui8matrix_erosion_SIMD_RR_row                    , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SIMD},
                                        {.func_name = "ui8matrix_erosion_erosion_SIMD_FO"             , .vec_morpho_func = ui8matrix_erosion_erosion_SIMD_FO                 , .pack_type = NO_PACK, .op_type = FUSION, .instr_type = SIMD},
                                        {.func_name = "ui8matrix_erosion_erosion_SIMD_FO_RR_row"      , .vec_morpho_func = ui8matrix_erosion_erosion_SIMD_FO_RR_row          , .pack_type = NO_PACK, .op_type = FUSION, .instr_type = SIMD},
                                        {.func_name = "ui8matrix_erosion_naive"                   , .morpho_func = ui8matrix_erosion_naive                                    , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR,},
                                        {.func_name = "ui8matrix_erosion_LU3x3_O1xO1"                , .morpho_func = ui8matrix_erosion_LU3x3_O1xO1                           , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR,},
                                        {.func_name = "ui8matrix_erosion_LU3x3_ExLU_O3"                , .morpho_func = ui8matrix_erosion_LU3x3_ExLU_O3                       , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR,},
                                        {.func_name = "ui8matrix_erosion_LU3x3_InLU_O3"                , .morpho_func = ui8matrix_erosion_LU3x3_InLU_O3                       , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR,},
                                        {.func_name = "ui8matrix_erosion_LU3x3_ComLU_O3"                , .morpho_func = ui8matrix_erosion_LU3x3_ComLU_O3                     , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR,},
                                        {.func_name = "ui8matrix_erosion_LU3x3_ExLU_O3_ValAddrRR"      , .morpho_func = ui8matrix_erosion_LU3x3_ExLU_O3_ValAddrRR             , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR,},
                                        {.func_name = "ui8matrix_erosion_LU3x3_InLU_O3_ValAddrRR"      , .morpho_func = ui8matrix_erosion_LU3x3_InLU_O3_ValAddrRR             , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR,},
                                        {.func_name = "ui8matrix_erosion_LU3x3_ComLU_O3_ValAddrRR"      , .morpho_func = ui8matrix_erosion_LU3x3_ComLU_O3_ValAddrRR           , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR,},
                                        {.func_name = "ui8matrix_erosion_LU3x3_InLU_O3_AddrRR"         , .morpho_func = ui8matrix_erosion_LU3x3_InLU_O3_AddrRR                , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR,},
                                        {.func_name = "ui8matrix_erosion_LU3x3_ExLU_O3_AddrRR"         , .morpho_func = ui8matrix_erosion_LU3x3_ExLU_O3_AddrRR                , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR,},
                                        {.func_name = "ui8matrix_erosion_LU3x3_ComLU_O3_AddrRR"         , .morpho_func = ui8matrix_erosion_LU3x3_ComLU_O3_AddrRR              , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR,},
                                        {.func_name = "ui8matrix_erosion_LU3x3_ExLU_O3_NS"             , .morpho_func = ui8matrix_erosion_LU3x3_ExLU_O3_NS                    , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR,},
                                        {.func_name = "ui8matrix_erosion_LU3x3_ExLU_O3_RR_NS"          , .morpho_func = ui8matrix_erosion_LU3x3_ExLU_O3_RR_NS                 , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR,},
                                        {.func_name = "ui8matrix_erosion_LU3x3_InLU_O3_NS"             , .morpho_func = ui8matrix_erosion_LU3x3_InLU_O3_NS                    , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR,},
                                        {.func_name = "ui8matrix_erosion_LU3x3_InLU_O3_RR_NS"          , .morpho_func = ui8matrix_erosion_LU3x3_InLU_O3_RR_NS                 , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR,},
                                        {.func_name = "ui8matrix_erosion_LU3x3_ComLU_O3_NS"             , .morpho_func = ui8matrix_erosion_LU3x3_ComLU_O3_NS                  , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR,},
                                        {.func_name = "ui8matrix_erosion_row_pipeline"               , .morpho_func = ui8matrix_erosion_row_pipeline                          , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR,},
                                        {.func_name = "ui8matrix_erosion_col_pipeline"               , .morpho_func = ui8matrix_erosion_col_pipeline                          , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR,},
                                        {.func_name = "ui8matrix_erosion_col_pipeline_RR"            , .morpho_func = ui8matrix_erosion_col_pipeline_RR                       , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR,},
                                        {.func_name = "ui8matrix_erosion_pipeline_LU3x3_ExLU_O3"       , .morpho_func = ui8matrix_erosion_pipeline_LU3x3_ExLU_O3              , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR,},
                                        {.func_name = "ui8matrix_erosion_pipeline_LU3x3_ExLU_O3_RR"    , .morpho_func = ui8matrix_erosion_pipeline_LU3x3_ExLU_O3_RR           , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR,},
                                        {.func_name = "ui8matrix_erosion_pipeline_LU3x3_InLU_O3"       , .morpho_func = ui8matrix_erosion_pipeline_LU3x3_InLU_O3              , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR,},
                                        {.func_name = "ui8matrix_erosion_pipeline_LU3x3_InLU_O3_RR"    , .morpho_func = ui8matrix_erosion_pipeline_LU3x3_InLU_O3_RR           , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR,},
                                        {.func_name = "ui8matrix_erosion_pipeline2_LU3x3_ExLU_O3"             , .morpho_func = ui8matrix_erosion_pipeline2_LU3x3_ExLU_O3      , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR,},
                                        {.func_name = "ui8matrix_erosion_pipeline2_LU3x3_ExLU_O3_RR"          , .morpho_func = ui8matrix_erosion_pipeline2_LU3x3_ExLU_O3_RR   , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR,},
                                        {.func_name = "ui8matrix_erosion_pipeline2_LU3x3_InLU_O3"             , .morpho_func = ui8matrix_erosion_pipeline2_LU3x3_InLU_O3      , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR,},
                                        {.func_name = "ui8matrix_erosion_pipeline2_LU3x3_InLU_O3_RR"          , .morpho_func = ui8matrix_erosion_pipeline2_LU3x3_InLU_O3_RR   , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR,},
                                        {.func_name = "ui8matrix_erosion_row_so"            , .morpho_func = ui8matrix_erosion_row_so         , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR,},
                                        {.func_name = "ui8matrix_erosion_row_so_InLU_O3"    , .morpho_func = ui8matrix_erosion_row_so_InLU_O3 , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR,},
                                        {.func_name = "ui8matrix_erosion_row_so_ExLU_O3"    , .morpho_func = ui8matrix_erosion_row_so_ExLU_O3 , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR,},
                                        {.func_name = "ui8matrix_erosion_col_so"            , .morpho_func = ui8matrix_erosion_col_so         , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR,},
                                        {.func_name = "ui8matrix_erosion_col_so_InLU_O3"    , .morpho_func = ui8matrix_erosion_col_so_InLU_O3 , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR,},
                                        {.func_name = "ui8matrix_erosion_col_so_ExLU_O3"    , .morpho_func = ui8matrix_erosion_col_so_ExLU_O3          , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR,},
                                        {.func_name = "ui8matrix_erosion_LU3x3_O1xO1_OMP"            , .morpho_func = ui8matrix_erosion_LU3x3_O1xO1_OMP                                , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR,},
                                        {.func_name = "ui8matrix_erosion_LU3x3_ExLU_O3_OMP"            , .morpho_func = ui8matrix_erosion_LU3x3_ExLU_O3_OMP                            , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR,},
                                        {.func_name = "ui8matrix_erosion_LU3x3_InLU_O3_OMP"            , .morpho_func = ui8matrix_erosion_LU3x3_InLU_O3_OMP                            , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR,},
                                        {.func_name = "ui8matrix_erosion_LU3x3_ComLU_O3_OMP"            , .morpho_func = ui8matrix_erosion_LU3x3_ComLU_O3_OMP                          , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR,},
                                        {.func_name = "ui8matrix_erosion_LU3x3_InLU_O3_AddrRR_OMP"     , .morpho_func = ui8matrix_erosion_LU3x3_InLU_O3_AddrRR_OMP                     , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR,},
                                        {.func_name = "ui8matrix_erosion_LU3x3_ExLU_O3_AddrRR_OMP"     , .morpho_func = ui8matrix_erosion_LU3x3_ExLU_O3_AddrRR_OMP                     , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR,},
                                        {.func_name = "ui8matrix_erosion_LU3x3_ComLU_O3_AddrRR_OMP"     , .morpho_func = ui8matrix_erosion_LU3x3_ComLU_O3_AddrRR_OMP                   , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR,},
                                        {.func_name = "ui8matrix_erosion_pipeline_LU3x3_ExLU_O3_RR_OMP"         , .morpho_func = ui8matrix_erosion_pipeline_LU3x3_ExLU_O3_RR_OMP       , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR,},
                                       {.func_name = "ui8matrix_erosion_pipeline2_LU3x3_InLU_O3_RR_OMP"      , .morpho_func = ui8matrix_erosion_pipeline2_LU3x3_InLU_O3_RR_OMP        , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR,},
                                        {.func_name = "ui8matrix_erosion_row_so_OMP"            , .morpho_func = ui8matrix_erosion_row_so_OMP          , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR,},
                                        {.func_name = "ui8matrix_erosion_row_so_InLU_O3_OMP"    , .morpho_func = ui8matrix_erosion_row_so_InLU_O3_OMP  , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR,},
                                        {.func_name = "ui8matrix_erosion_row_so_ExLU_O3_OMP"    , .morpho_func = ui8matrix_erosion_row_so_ExLU_O3_OMP  , .pack_type = NO_PACK, .op_type = NORMAL, .instr_type = SCALAR,},
                                       };
    
    struct sd_set SDs_step0[] = {//SigmaDelta_step0_InLU_O3_OMP
                                {.func_name = "SigmaDelta_step0_mem", .sd_func = SigmaDelta_step0_mem, .n_coeff = N, .v_min = Vmin, .v_max=Vmax, .instr_type = SCALAR},
                                {.func_name = "SigmaDelta_step0_naive", .sd_func = SigmaDelta_step0_naive, .n_coeff = N, .v_min = Vmin, .v_max=Vmax, .instr_type = SCALAR},
                                {.func_name = "SigmaDelta_step0_SIMD", .vec_sd_func = SigmaDelta_step0_SIMD, .n_coeff = N, .v_min = Vmin, .v_max=Vmax, .instr_type = SIMD},
                                {.func_name = "SigmaDelta_step0_SIMD_and_load", .vec_sd_func = SigmaDelta_step0_SIMD_and_load, .n_coeff = N, .v_min = Vmin, .v_max=Vmax, .instr_type = SIMD},
                                {.func_name = "SigmaDelta_step0_SIMD_or_load", .vec_sd_func = SigmaDelta_step0_SIMD_or_load, .n_coeff = N, .v_min = Vmin, .v_max=Vmax, .instr_type = SIMD},
                                {.func_name = "SigmaDelta_step0_SIMD_load_load", .vec_sd_func = SigmaDelta_step0_SIMD_load_load, .n_coeff = N, .v_min = Vmin, .v_max=Vmax, .instr_type = SIMD},
                                {.func_name = "SigmaDelta_step0_SIMD_store_load", .vec_sd_func = SigmaDelta_step0_SIMD_store_load, .n_coeff = N, .v_min = Vmin, .v_max=Vmax, .instr_type = SIMD},
                                {.func_name = "SigmaDelta_step0_SIMD_memset_load", .vec_sd_func = SigmaDelta_step0_SIMD_memset_load, .n_coeff = N, .v_min = Vmin, .v_max=Vmax, .instr_type = SIMD},
                                {.func_name = "SigmaDelta_step0_InLU_O3_OMP", .vec_sd_func = SigmaDelta_step0_SIMD_InLU_O3_OMP, .n_coeff = N, .v_min = Vmin, .v_max=Vmax, .instr_type = SIMD},
                                {.func_name = "SigmaDelta_step0_InLU_O4_OMP", .vec_sd_func = SigmaDelta_step0_SIMD_InLU_O4_OMP, .n_coeff = N, .v_min = Vmin, .v_max=Vmax, .instr_type = SIMD},
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
                                {.func_name = "SigmaDelta_step2_InLU_O3_SIMD", .vec_sd_func = SigmaDelta_step2_InLU_O3_SIMD, .n_coeff = N, .v_min = Vmin, .v_max=Vmax, .instr_type = SIMD},
                                {.func_name = "SigmaDelta_step2_ExLU_O3_SIMD", .vec_sd_func = SigmaDelta_step2_ExLU_O3_SIMD, .n_coeff = N, .v_min = Vmin, .v_max=Vmax, .instr_type = SIMD},
                                }
                                ;
    struct sd_set SDs_step3[] = {
                                {.func_name = "SigmaDelta_step3_naive", .sd_func = SigmaDelta_step3_naive, .n_coeff = N, .v_min = Vmin, .v_max =Vmax, .instr_type = SCALAR},
                                {/*Best :2.94*/.func_name = "SigmaDelta_step3_InLU_O3_NoIf", .sd_func = SigmaDelta_step3_InLU_O3_NoIf, .n_coeff = N, .v_min = Vmin, .v_max =Vmax, .instr_type = SCALAR},
                                {              .func_name = "SigmaDelta_step3_ExLU_O3_NoIf", .sd_func = SigmaDelta_step3_ExLU_O3_NoIf, .n_coeff = N, .v_min = Vmin, .v_max =Vmax, .instr_type = SCALAR},
                                {.func_name = "SigmaDelta_step3_SIMD", .vec_sd_func = SigmaDelta_step3_SIMD, .n_coeff = N, .v_min = Vmin, .v_max=Vmax, .instr_type = SIMD},
                                {.func_name = "SigmaDelta_step3_SIMD_ver2", .vec_sd_func = SigmaDelta_step3_SIMD_ver2, .n_coeff = N, .v_min = Vmin, .v_max=Vmax, .instr_type = SIMD},
                                {.func_name = "SigmaDelta_step3_InLU_O3_SIMD", .vec_sd_func = SigmaDelta_step3_InLU_O3_SIMD, .n_coeff = N, .v_min = Vmin, .v_max=Vmax, .instr_type = SIMD},
                                };//
    struct sd_set SDs_step4[] = {
                                {.func_name = "SigmaDelta_step4_naive", .sd_func = SigmaDelta_step4_naive, .n_coeff = N, .v_min = Vmin, .v_max =Vmax, .instr_type = SCALAR},
                                {/*Best :0.18*/.func_name = "SigmaDelta_step4_InLU_O3_NoIf", .sd_func = SigmaDelta_step4_InLU_O3_NoIf, .n_coeff = N, .v_min = Vmin, .v_max =Vmax, .instr_type = SCALAR},
                                {              .func_name = "SigmaDelta_step4_ExLU_O3_NoIf", .sd_func = SigmaDelta_step4_ExLU_O3_NoIf, .n_coeff = N, .v_min = Vmin, .v_max =Vmax, .instr_type = SCALAR},
                                {.func_name = "SigmaDelta_step4_SIMD", .vec_sd_func = SigmaDelta_step4_SIMD, .n_coeff = N, .v_min = Vmin, .v_max=Vmax, .instr_type = SIMD},
                                };
                                
                                
    struct sd_set SD_steps[] =  {
                                {.func_name = "SigmaDelta_step0_naive", .sd_func = SigmaDelta_step0_naive, .n_coeff = N, .v_min = Vmin, .v_max=Vmax, .instr_type = SCALAR},
                                {.func_name = "SigmaDelta_step0_mem", .sd_func = SigmaDelta_step0_mem, .n_coeff = N, .v_min = Vmin, .v_max=Vmax, .instr_type = SCALAR},
                                {.func_name = "SigmaDelta_step0_SIMD_memset_load", .vec_sd_func = SigmaDelta_step0_SIMD_memset_load, .n_coeff = N, .v_min = Vmin, .v_max=Vmax, .instr_type = SIMD},
                                {.func_name = "SigmaDelta_step1_naive"  , .sd_func = SigmaDelta_step1_naive, .n_coeff = N, .v_min = Vmin, .v_max=Vmax, .instr_type = SCALAR},
                                {.func_name = "SigmaDelta_step1_InLU_O3_NoIf", .sd_func = SigmaDelta_step1_InLU_O3_NoIf, .n_coeff = N, .v_min = Vmin, .v_max=Vmax, .instr_type = SCALAR},
                                {.func_name = "SigmaDelta_step1_ExLU_O3_NoIf", .sd_func = SigmaDelta_step1_ExLU_O3_NoIf, .n_coeff = N, .v_min = Vmin, .v_max=Vmax, .instr_type = SCALAR},
                                {.func_name = "SigmaDelta_step1_SIMD", .vec_sd_func = SigmaDelta_step1_SIMD, .n_coeff = N, .v_min = Vmin, .v_max=Vmax, .instr_type = SIMD},
                                {.func_name = "SigmaDelta_step2_naive", .sd_func = SigmaDelta_step2_naive, .n_coeff = N, .v_min = Vmin, .v_max =Vmax, .instr_type = SCALAR},
                                {.func_name = "SigmaDelta_step2_InLU_O3_bitop", .sd_func = SigmaDelta_step2_InLU_O3_bitop, .n_coeff = N, .v_min = Vmin, .v_max=Vmax, .instr_type = SCALAR},
                                {.func_name = "SigmaDelta_step2_ExLU_O3_bitop", .sd_func = SigmaDelta_step2_ExLU_O3_bitop, .n_coeff = N, .v_min = Vmin, .v_max=Vmax, .instr_type = SCALAR},
                                {.func_name = "SigmaDelta_step2_SIMD", .vec_sd_func = SigmaDelta_step2_SIMD, .n_coeff = N, .v_min = Vmin, .v_max=Vmax, .instr_type = SIMD},
                                {.func_name = "SigmaDelta_step2_InLU_O3_SIMD", .vec_sd_func = SigmaDelta_step2_InLU_O3_SIMD, .n_coeff = N, .v_min = Vmin, .v_max=Vmax, .instr_type = SIMD},
                                {.func_name = "SigmaDelta_step2_ExLU_O3_SIMD", .vec_sd_func = SigmaDelta_step2_ExLU_O3_SIMD, .n_coeff = N, .v_min = Vmin, .v_max=Vmax, .instr_type = SIMD},
                                {.func_name = "SigmaDelta_step3_naive", .sd_func = SigmaDelta_step3_naive, .n_coeff = N, .v_min = Vmin, .v_max =Vmax, .instr_type = SCALAR},
                                {.func_name = "SigmaDelta_step3_InLU_O3_NoIf", .sd_func = SigmaDelta_step3_InLU_O3_NoIf, .n_coeff = N, .v_min = Vmin, .v_max =Vmax, .instr_type = SCALAR},
                                {.func_name = "SigmaDelta_step3_ExLU_O3_NoIf", .sd_func = SigmaDelta_step3_ExLU_O3_NoIf, .n_coeff = N, .v_min = Vmin, .v_max =Vmax, .instr_type = SCALAR},
                                {.func_name = "SigmaDelta_step3_SIMD", .vec_sd_func = SigmaDelta_step3_SIMD, .n_coeff = N, .v_min = Vmin, .v_max=Vmax, .instr_type = SIMD},
                                {.func_name = "SigmaDelta_step3_InLU_O3_SIMD", .vec_sd_func = SigmaDelta_step3_InLU_O3_SIMD, .n_coeff = N, .v_min = Vmin, .v_max=Vmax, .instr_type = SIMD},
                                {.func_name = "SigmaDelta_step4_naive", .sd_func = SigmaDelta_step4_naive, .n_coeff = N, .v_min = Vmin, .v_max =Vmax, .instr_type = SCALAR},
                                {.func_name = "SigmaDelta_step4_InLU_O3_NoIf", .sd_func = SigmaDelta_step4_InLU_O3_NoIf, .n_coeff = N, .v_min = Vmin, .v_max =Vmax, .instr_type = SCALAR},
                                {.func_name = "SigmaDelta_step4_ExLU_O3_NoIf", .sd_func = SigmaDelta_step4_ExLU_O3_NoIf, .n_coeff = N, .v_min = Vmin, .v_max =Vmax, .instr_type = SCALAR},
                                {.func_name = "SigmaDelta_step4_SIMD", .vec_sd_func = SigmaDelta_step4_SIMD, .n_coeff = N, .v_min = Vmin, .v_max=Vmax, .instr_type = SIMD},
                                };
    
    struct complete_sd_set completeSDs[] = {//
                                    {.func_name = "SigmaDelta_naive", .sd_step0 = SigmaDelta_step0_naive, .sd_func = SigmaDelta_naive, .n_coeff = N, .v_min = Vmin, .v_max = Vmax, .instr_type = SCALAR},                                  
                                    {.func_name = "SigmaDelta_best", .sd_step0 = SigmaDelta_step0_mem, .sd_func = SigmaDelta_best, .n_coeff = N, .v_min = Vmin, .v_max = Vmax, .instr_type = SCALAR},
                                    {.func_name = "SigmaDelta_SIMD", .vec_sd_step0 = SigmaDelta_step0_SIMD, .vec_sd_func = SigmaDelta_SIMD, .n_coeff = N, .v_min = Vmin, .v_max = Vmax, .instr_type = SIMD},
                                    {.func_name = "SigmaDelta_SIMD_FL", .vec_sd_step0 = SigmaDelta_step0_SIMD_memset_load, .vec_sd_func = SigmaDelta_SIMD_FL, .n_coeff = N, .v_min = Vmin, .v_max = Vmax, .instr_type = SIMD},
                                    {.func_name = "SigmaDelta_SIMD_FL_OMP", .vec_sd_step0 = SigmaDelta_step0_SIMD_memset_load, .vec_sd_func = SigmaDelta_SIMD_FL_OMP, .n_coeff = N, .v_min = Vmin, .v_max = Vmax, .instr_type = SIMD},
                                  };
    struct morpho_set sequences[] = {
                                        //{.func_name = "ui8matrix_sequence_SIMD_Pipeline_FO_InLU_O3_ValAddrRR_OMP"   , .vec_morpho_func = ui8matrix_sequence_SIMD_Pipeline_FO_InLU_O3_ValAddrRR_OMP     , .pack_type=NO_PACK , .op_type=NORMAL      , .instr_type = SIMD}, 
//                                        {.func_name = "ui8matrix_sequence_SIMD_Pipeline_FO_InLU_O3_ValAddrRR"       , .vec_morpho_func = ui8matrix_sequence_SIMD_Pipeline_FO_InLU_O3_ValAddrRR     , .pack_type=NO_PACK , .op_type=NORMAL      , .instr_type = SIMD}, 
                                        {.func_name = "ui8matrix_sequence_SIMD_Pipeline_FO_InLU_O3_ValAddrRR"       , .vec_morpho_func = ui8matrix_sequence_SIMD_Pipeline_FO_InLU_O3_ValAddrRR     , .pack_type=NO_PACK , .op_type=NORMAL      , .instr_type = SIMD}, 
                                        {.func_name = "ui8matrix_sequence_SIMD_FO_InLU_O3_ValAddrRR_OMP"            , .vec_morpho_func = ui8matrix_sequence_SIMD_FO_InLU_O3_ValAddrRR_OMP     , .pack_type=NO_PACK , .op_type=NORMAL      , .instr_type = SIMD}, 
                                        {.func_name = "ui8matrix_sequence_SIMD_FO_InLU_O3_ValAddrRR"                , .vec_morpho_func = ui8matrix_sequence_SIMD_FO_InLU_O3_ValAddrRR     , .pack_type=NO_PACK , .op_type=NORMAL      , .instr_type = SIMD}, 
                                        {.func_name = "ui8matrix_sequence_row_so_fo"                                  , .morpho_func = ui8matrix_sequence_row_so_fo, .pack_type=NO_PACK , .op_type=NORMAL      , .instr_type = SCALAR}, 
                                        {.func_name = "ui8matrix_sequence_naive"                                    , .morpho_func = ui8matrix_sequence_naive     , .pack_type=NO_PACK , .op_type=NORMAL      , .instr_type = SCALAR}, 
                                        {.func_name = "ui8matrix_sequence_col_so"                                     , .morpho_func = ui8matrix_sequence_col_so, .pack_type=NO_PACK , .op_type=NORMAL      , .instr_type = SCALAR}, 
                                        {.func_name = "ui8matrix_sequence_row_so"                                     , .morpho_func = ui8matrix_sequence_row_so, .pack_type=NO_PACK , .op_type=NORMAL      , .instr_type = SCALAR}, 
                                        {.func_name = "ui8matrix_sequence_row_so_fo_pipeline"                         , .morpho_func = ui8matrix_sequence_row_so_fo_pipeline, .pack_type=NO_PACK , .op_type=NORMAL      , .instr_type = SCALAR}, 
                                        {.func_name = "ui8matrix_sequence_row_so_fo_pipeline2"                        , .morpho_func = ui8matrix_sequence_row_so_fo_pipeline2, .pack_type=NO_PACK , .op_type=NORMAL      , .instr_type = SCALAR}, 
                                        // {.func_name = "ui8matrix_sequence_row_so_OMP"  , .morpho_func = ui8matrix_sequence_row_so_OMP   , .instr_type = SCALAR},
                                       };
    struct complete_process_set cps[] = {
                                        {.func_name = "SD_Naive/Sequential_Morpho_Naive"                 , .sd_step0 = SigmaDelta_step0_naive, .sd_func = SigmaDelta_naive, .morpho_func = ui8matrix_sequence_naive, .instr_type= SCALAR},
                                        {.func_name = "SD_FO+SIMD/Morpho-Pipeline+FO+InLU_O3+FullRR+SIMD", .vec_sd_step0 = SigmaDelta_step0_SIMD_memset_load, .vec_sd_func =SigmaDelta_SIMD_FL, .vec_morpho_func = ui8matrix_sequence_SIMD_Pipeline_FO_InLU_O3_ValAddrRR, .instr_type= SIMD},
                                        {.func_name = "SD_FO+SIMD+OMP/Morpho-Pipeline+FO+InLU_O3+FullRR+SIMD ", .vec_sd_step0 = SigmaDelta_step0_SIMD_memset_load, .vec_sd_func =SigmaDelta_SIMD_FL_OMP, .vec_morpho_func = ui8matrix_sequence_SIMD_Pipeline_FO_InLU_O3_ValAddrRR, .instr_type= SIMD},
                                        {.func_name = "SD_FO+SIMD+OMP/Morpho-FO+InLU_O3+FullRR+SIMD+OMP ", .vec_sd_step0 = SigmaDelta_step0_SIMD_memset_load, .vec_sd_func =SigmaDelta_SIMD_FL_OMP, .vec_morpho_func = ui8matrix_sequence_SIMD_FO_InLU_O3_ValAddrRR_OMP, .instr_type= SIMD},
                                        };
void launch_movement_detection(char *filename_format, int start, int end, char *res_filename_format)
{
    vuint8 **vTempBuffer;
    uint8 **ui8image;
    p_vimage t0, t1;
    char filename[128], res_filename[128];
    long nrl, nrh, ncl, nch;
    int cnt = 0, error_no = 0;


    sprintf(filename, filename_format, start);
    t0 = create_vimage(filename);
    
    SigmaDelta_step0_SIMD_memset_load(t0->M, t0->I, t0->V, t0->nrl, t0->nrh, t0->v0, t0->v1, N, Vmin, Vmax);
    vTempBuffer = vui8matrix(t0->nrl, t0->nrh, t0->v0, t0->v1);
    
    for (int i = start + 1; i < end + 1; i++, cnt++) {
        sprintf(filename, filename_format, i);
        sprintf(res_filename, res_filename_format, cnt);
        t1 = create_vimage(filename);
        SigmaDelta_SIMD_FL(t0, t1, N, Vmin, Vmax);
        ui8matrix_sequence_SIMD_Pipeline_FO_InLU_O3_ValAddrRR(t1->E, t1->nrl + BORD, t1->nrh - BORD, t1->ncl + BORD, t1->nch - BORD, t1->v0 + vBORD, t1->v1 -vBORD, vTempBuffer, t1->Omega);
        ui8image = vui8matrix_to_ui8matrix(t1->Omega, t1->nrl + BORD, t1->nrh - BORD, t1->v0 + vBORD, t1->v1 - vBORD, &nrl, &nrh, &ncl, &nch);
        binary_to_octal_ui8matrix(ui8image, t1->nrl + BORD, t1->nrh - BORD, t1->ncl + BORD, t1->nch + BORD);
        SavePGM_ui8matrix(ui8image, t1->nrl + BORD, t1->nrh - BORD, t1->ncl + BORD, t1->nch + BORD,  res_filename);
        free_vimage(t0);
        free_ui8matrix(ui8image, nrl, nrh, ncl, nch);
        t0 = t1;
    }
    free_vui8matrix(vTempBuffer, t0->nrl, t0->nrh, t0->v0, t0->v1);
}

void launch_naive_movement_detection(char *filename_format, int start, int end, char *res_filename_format)
{
    uint8 **tempBuffer;
    p_image t0, t1;
    char filename[128], res_filename[128];
    long nrl, nrh, ncl, nch;
    int cnt = 0, error_no = 0;


    sprintf(filename, filename_format, start);
    t0 = create_image(filename);
    
    SigmaDelta_step0_naive(t0->M, t0->I, t0->V, t0->nrl, t0->nrh, t0->ncl, t0->nch, N, Vmin, Vmax);
   
    tempBuffer = ui8matrix(t0->nrl, t0->nrh, t0->ncl, t0->nch);
    for (int i = start + 1; i < end + 1; i++, cnt++) {
        sprintf(filename, filename_format, i);
        sprintf(res_filename, res_filename_format, cnt);
        t1 = create_image(filename);
        SigmaDelta_naive(t0, t1, N, Vmin, Vmax); 
        ui8matrix_sequence_naive(t1->E, t1->nrl + BORD, t1->nrh - BORD, t1->ncl + BORD, t1->nch - BORD, tempBuffer, t1->Omega);
        binary_to_octal_ui8matrix(t1->Omega, t1->nrl + BORD, t1->nrh - BORD, t1->ncl + BORD, t1->nch + BORD);
        SavePGM_ui8matrix(t1->Omega, t1->nrl, t1->nrh, t1->ncl, t1->nch,  res_filename);
        free_image(t0);
        t0 = t1;
    }
    free_ui8matrix(tempBuffer, t0->nrl, t0->nrh, t0->ncl, t0->nch);
    free_image(t1);
}
// -----------
int main(void)
// -----------
{
    
    long nb_sets;

    
    // launch_movement_detection("../car3/car_%d.pgm", 3000, 3199, "../res/res%04d_a.pgm");
    // launch_naive_movement_detection("../car3/car_%d.pgm", 3000, 3199, "../res/res%04d_b.pgm");
    //test_SigmaDelta_step0("../car3/car_3000.pgm", "../car3/car_3001.pgm", SDs_step0, 10, false);
    //test_SigmaDelta_step1("../car3/car_3000.pgm", "../car3/car_3001.pgm", SDs_step1, 4, false);
    //test_SigmaDelta_step2("../car3/car_3000.pgm", "../car3/car_3001.pgm", SDs_step2, 6, false);
    //test_SigmaDelta_step3("../car3/car_3000.pgm", "../car3/car_3001.pgm", SDs_step3, 6, false);
    //test_SigmaDelta_step4("../car3/car_3000.pgm", "../car3/car_3001.pgm", SDs_step4, 4, false);
    //test_SigmaDelta("../car3/car_3000.pgm", "../car3/car_3001.pgm", completeSDs, 4, false);
    
    // launch_SD_step_benchmark("output/benchmark_sdstep0.dat"             , SDs_step0, 10, 1, 1, 100, 10000, 100);
    // launch_SD_step_benchmark("output/benchmark_sdstep1.dat"             , SDs_step1, 4, 1, 1, 200, 5000, 100);
    // launch_SD_step_benchmark("output/benchmark_sdstep2.dat"             , SDs_step2, 6, 1, 1, 200, 5000, 50);
    // launch_SD_step_benchmark("output/benchmark_sdstep3.dat"             , SDs_step3, 5, 1, 1, 200, 5000, 50);
    // launch_SD_step_benchmark("output/benchmark_sdstep4.dat"             , SDs_step4, 4, 1, 1, 200, 5000, 100);
    // launch_SD_step_benchmark("output/benchmark_SD_step.dat"       , SD_steps   ,       5, 1, 1, 200, 10000, 100);
    // launch_SD_benchmark(     "output/benchmark_SD.dat"            , completeSDs, 1, 1, 1, 100, 10000, 100);
    


    //test_dilations("../car3/car_3000.pgm", dilations, 53, false);
//    test_erosions ("../car3/car_3000.pgm", erosions , 56, false);
    //test_sequences ("../car3/car_3001.pgm", sequences , 10, false);
  //  launch_morpho_benchmark( "output/benchmark_dilation.csv", dilations  , 51, 1, 1,  100, 32000, 32);
//    launch_morpho_benchmark( "output/benchmark_erosion.csv" , erosions   , 54, 1, 1,  100, 32000, 32);
    launch_morpho_benchmark( "output/benchmark_sequence.csv" , sequences , 9, 1, 1,  32, 6100, 32);
    launch_complete_process_benchmark("output/full_benchmark.csv", cps   ,  4, 1, 1, 32, 6100, 32);
    


    
    
    return 0;    
}

























