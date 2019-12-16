#ifndef  __TEST_MORPHO_H__
#define __TEST_MORPHO_H__

#pragma message("  include  test_morpho.h")

// A struct type that contains a morpho function and a structuring element
struct morpho_set{
    char func_name[128];
    void (*morpho_func)(uint8** ppInput, long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput);
    struct struct_elem_dim *s;
};


double **benchmark(struct morpho_set *msets, long nb_sets, long ls, long hs, long step, int nb_tests);
unsigned long long get_cpu_cycles(struct morpho_set *ptr_mset, uint8 **ppInput,  long nrl, long nrh, long ncl, long nch, uint8 **ppOutput);
unsigned long long get_min_cpu_cycles(struct morpho_set *ptr_mset, long packet_size, uint8 **ppInput,  long nrl, long nrh, long ncl, long nch, uint8 **ppOutput);


uint8** prologue_test_integration_3x3(struct morpho_set *mset, uint8** ppInput, uint8** ppOutput);
uint8** prologue_test_integration_5x5(struct morpho_set *mset, uint8** ppInput, uint8** ppOutput);
void epilogue_test_integration_3x3(struct morpho_set *mset, uint8** ppInput, uint8** ppOutput);
void epilogue_test_integration_5x5(struct morpho_set *mset, uint8** ppInput, uint8** ppOutput);

void test_implementation_erosion_3x3 (struct morpho_set *erosion_set );
void test_implementation_erosion_5x5 (struct morpho_set *erosion_set );
void test_implementation_dilation_3x3(struct morpho_set *dilation_set);
void test_implementation_dilation_5x5(struct morpho_set *dilation_set);

void test_integration_erosion_3x3 (struct morpho_set *erosion_set);
void test_integration_erosion_5x5 (struct morpho_set *erosion_set);
void test_integration_dilation_3x3(struct morpho_set *dilation_set);
void test_integration_dilation_5x5(struct morpho_set *dilation_set);

bool morpho_produces_one(struct morpho_set *mset, uint8** ppInput);

bool check_dimension_of_square_structuring_element(p_struct_elem_dim s,  long size);

static inline bool check_for_3x3_structuring_element(p_struct_elem_dim s);
static inline bool check_for_3x3_structuring_element(p_struct_elem_dim s) 
{
    return check_dimension_of_square_structuring_element(s, 3);
}
static inline bool check_for_5x5_structuring_element(p_struct_elem_dim s);
static inline bool check_for_5x5_structuring_element(p_struct_elem_dim s) 
{
    return check_dimension_of_square_structuring_element(s, 5);
}


static inline uint32 extract_bits_from_permutation(uint32 perm, long irow, long ncol);
static inline uint32 extract_bits_from_permutation(uint32 perm, long irow, long ncol) {
    return (perm >> irow * ncol);
}

static inline uint8 get_column_at(uint8 col_vals, long icol);
static inline uint8 get_column_at(uint8 col_vals, long icol) {
    return (col_vals >> icol) & 0x1;
}

uint8 **rand_ui8matrix(long size, p_struct_elem_dim s);
uint8 **ui8matrix_checker(long nrl, long nrh, long ncl, long nch, const long chkr_size, const uint8 xor_mask);

double **init_benchmark_results(long nb_funcs, long size);
void     free_benchmark_results(double **results, long nb_funcs);

uint8**  ui8matrix_permutation (uint8** m, long nrl, long nrh, long ncl, long nch, uint32 perm);


#endif /*  __TEST_MORPHO_H__ */