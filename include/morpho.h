/* --------------- */
/* --- morpho.h --- */
/* --------------- */
#ifndef __MORPHO_H__
#define __MORPHO_H__

#pragma message("  include  morpho.h")

// A struct type that contains a morpho function and a structuring element
struct morpho_set{
    char func_name[128];
    void (*morpho_func)(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
    void (*vec_morpho_func)(vuint8** vX, int i0, int i1, int j0, int j1, vuint8 **vTempBuffer, vuint8 **vY);
    enum {NO_PACK, HPACK, VPACK}pack_type;
    enum {NORMAL, FUSION}op_type;
    instruction_type instr_type;
};



// typedef void (*morpho_func_t)(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);

void ui8matrix_dilation_naive(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_dilation5_naive(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
/****************************/
/******* Loop Unroll ********/
/****************************/
void ui8matrix_dilation_LU3x3_O1xO1 (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_dilation_LU3x3_ExLU_O3 (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_dilation_LU3x3_InLU_O3 (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_dilation_LU3x3_ComLU_O3 (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);

/***********************************************************************/
/******** Loop Unroll + Register Rotation of Values / Addresses ********/
/***********************************************************************/
void ui8matrix_dilation_LU3x3_ExLU_O3_ValAddrRR (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_dilation_LU3x3_InLU_O3_ValAddrRR (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_dilation_LU3x3_ComLU_O3_ValAddrRR (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);

/**************************************************************/
/******** Loop Unroll + Register Rotation of Addresses ********/
/**************************************************************/
void ui8matrix_dilation_LU3x3_InLU_O3_AddrRR (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_dilation_LU3x3_ExLU_O3_AddrRR (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_dilation_LU3x3_ComLU_O3_AddrRR (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);

void ui8matrix_dilation_pipeline2_LU3x3_InLU_O3 (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_dilation_pipeline2_LU3x3_ExLU_O3 (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_dilation_pipeline2_LU3x3_InLU_O3_RR (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_dilation_pipeline2_LU3x3_ExLU_O3_RR (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);


// void ui8matrix_dilation_LU3x3_InLU_O3_Full_Scalari (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
/************************************/
/******** Alternate Versions ********/
/************************************/
void ui8matrix_dilation5_divide_row_and_conquer (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **Y, uint8 **Z);
void ui8matrix_dilation_divide_row_and_conquer (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_dilation_divide_row_and_conquer_ExLU_O3 (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_dilation_divide_row_and_conquer_InLU_O3 (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_dilation_divide_col_and_conquer (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_dilation_divide_col_and_conquer_ExLU_O3 (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_dilation_divide_col_and_conquer_InLU_O3 (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);

void ui8matrix_dilation_LU3x3_ExLU_O3_NS (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_dilation_LU3x3_ExLU_O3_RR_NS (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_dilation_LU3x3_InLU_O3_NS (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_dilation_LU3x3_InLU_O3_RR_NS (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_dilation_LU3x3_ComLU_O3_RR_ver1 (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_dilation_LU3x3_ComLU_O3_NS (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);

/*************************/
/******* Pipeline ********/
/*************************/
void ui8matrix_dilation_row_pipeline (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_dilation_col_pipeline (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_dilation_col_pipeline_RR (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_dilation_pipeline_LU3x3_ExLU_O3 (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_dilation_pipeline_LU3x3_ExLU_O3_RR (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);

void ui8matrix_dilation_pipeline_LU3x3_InLU_O3 (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_dilation_pipeline_LU3x3_InLU_O3_RR (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);

/**************************/
/********* openMP *********/
/**************************/
void ui8matrix_dilation_LU3x3_O1xO1_OMP			   (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_dilation_LU3x3_ExLU_O3_OMP			   (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_dilation_LU3x3_InLU_O3_OMP			   (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_dilation_LU3x3_ComLU_O3_OMP			   (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_dilation_LU3x3_InLU_O3_AddrRR_OMP     (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_dilation_LU3x3_ExLU_O3_AddrRR_OMP     (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_dilation_LU3x3_ComLU_O3_AddrRR_OMP     (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_dilation_pipeline2_LU3x3_InLU_O3_OMP         (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_dilation_pipeline2_LU3x3_ExLU_O3_OMP         (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_dilation_pipeline2_LU3x3_InLU_O3_RR_OMP      (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_dilation_pipeline2_LU3x3_ExLU_O3_RR_OMP      (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_dilation_divide_row_and_conquer_OMP (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_dilation_divide_row_and_conquer_ExLU_O3_OMP (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_dilation_divide_row_and_conquer_InLU_O3_OMP (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_dilation_pipeline2_LU3x3_InLU_O3_RR_OMP (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_dilation_pipeline_LU3x3_ExLU_O3_RR_OMP (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);




void ui8matrix_erosion_naive(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_erosion5_naive(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
/****************************/
/******* Loop Unroll ********/
/****************************/
void ui8matrix_erosion_LU3x3_O1xO1 (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_erosion_LU3x3_ExLU_O3 (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_erosion_LU3x3_InLU_O3 (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_erosion_LU3x3_ComLU_O3 (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);

/***********************************************************************/
/******** Loop Unroll + Register Rotation of Values / Addresses ********/
/***********************************************************************/
void ui8matrix_erosion_LU3x3_ExLU_O3_ValAddrRR (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_erosion_LU3x3_InLU_O3_ValAddrRR (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_erosion_LU3x3_ComLU_O3_ValAddrRR (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);

/**************************************************************/
/******** Loop Unroll + Register Rotation of Addresses ********/
/**************************************************************/
void ui8matrix_erosion_LU3x3_InLU_O3_AddrRR (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_erosion_LU3x3_ExLU_O3_AddrRR (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_erosion_LU3x3_ComLU_O3_AddrRR (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);

void ui8matrix_erosion_pipeline2_LU3x3_InLU_O3 (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_erosion_pipeline2_LU3x3_ExLU_O3 (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_erosion_pipeline2_LU3x3_InLU_O3_RR (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_erosion_pipeline2_LU3x3_ExLU_O3_RR (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);


// void ui8matrix_erosion_LU3x3_InLU_O3_Full_Scalari (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
/************************************/
/******** Alternate Versions ********/
/************************************/

void ui8matrix_erosion_divide_row_and_conquer (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_erosion_divide_row_and_conquer_ExLU_O3 (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_erosion_divide_row_and_conquer_InLU_O3 (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_erosion_divide_col_and_conquer (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_erosion_divide_col_and_conquer_ExLU_O3 (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_erosion_divide_col_and_conquer_InLU_O3 (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);

void ui8matrix_erosion_LU3x3_ExLU_O3_NS (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_erosion_LU3x3_ExLU_O3_RR_NS (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_erosion_LU3x3_InLU_O3_NS (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_erosion_LU3x3_InLU_O3_RR_NS (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_erosion_LU3x3_ComLU_O3_RR_ver1 (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_erosion_LU3x3_ComLU_O3_NS (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);

/*************************/
/******* Pipeline ********/
/*************************/
void ui8matrix_erosion_row_pipeline (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_erosion_col_pipeline (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_erosion_col_pipeline_RR (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_erosion_pipeline_LU3x3_ExLU_O3 (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_erosion_pipeline_LU3x3_ExLU_O3_RR (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);

void ui8matrix_erosion_pipeline_LU3x3_InLU_O3 (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_erosion_pipeline_LU3x3_InLU_O3_RR (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);

/**************************/
/********* packed *********/
/**************************/
void ui8matrix_dilation_hpacked_divide_row_and_conquer (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_dilation_vpacked_divide_row_and_conquer (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);

/**************************/
/********* openMP *********/
/**************************/
void ui8matrix_erosion_LU3x3_O1xO1_OMP			         (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_erosion_LU3x3_ExLU_O3_OMP			     (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_erosion_LU3x3_InLU_O3_OMP			     (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_erosion_LU3x3_ComLU_O3_OMP			     (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_erosion_LU3x3_InLU_O3_AddrRR_OMP          (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_erosion_LU3x3_ExLU_O3_AddrRR_OMP          (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_erosion_LU3x3_ComLU_O3_AddrRR_OMP         (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_erosion_pipeline2_LU3x3_InLU_O3_OMP       (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_erosion_pipeline2_LU3x3_ExLU_O3_OMP       (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_erosion_pipeline2_LU3x3_InLU_O3_RR_OMP    (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_erosion_pipeline2_LU3x3_ExLU_O3_RR_OMP    (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_erosion_divide_row_and_conquer_OMP        (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_erosion_divide_row_and_conquer_ExLU_O3_OMP(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_erosion_divide_row_and_conquer_InLU_O3_OMP(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_erosion_pipeline2_LU3x3_InLU_O3_RR_OMP    (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_erosion_pipeline_LU3x3_ExLU_O3_RR_OMP     (uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);



void ui8matrix_sequence_naive(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_sequence_crnc(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **Y, uint8 **Z);
void ui8matrix_sequence_drnc(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **Y, uint8 **Z);
void ui8matrix_sequence_drnc_fo(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_sequence_LE_FO_pipeline(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **Y, uint8 **Z);
void ui8matrix_sequence_drnc_fo_pipeline(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_sequence_drnc_fo_pipeline2(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_sequence_fo(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
void ui8matrix_sequence_divide_row_and_conquer_OMP(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);




/******************/
/*** DEPRECATED ***/
/******************/

// uint8  erosion_naive(uint8** X, long row, long col);
// uint8 dilation_naive(uint8** X, long row, long col);



// void ui8matrix_lambda_morpho(morpho_func_t morpho, long order, uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
// void dilation_naive(uint8** X, long row, long col, uint8 **tempBuffer, uint8 **Y);
// void  erosion_naive(uint8** X, long row, long col, uint8 **tempBuffer, uint8 **Y);


void image_chain_processing(p_image img, int idx);
#endif /* __MORPHO_H__ */