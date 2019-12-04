#ifndef  __TEST_MORPHO_H__
#define __TEST_MORPHO_H__

#pragma message("  include  test_morpho.h")


uint8 **rand_ui8matrix(long size, p_struct_elem_dim s);
uint8 **ui8matrix_test_input(long size, const uint8 xor_mask, p_struct_elem_dim s);

double **init_benchmark_results(long nb_funcs, long size);
void     free_benchmark_results(double **results, long nb_funcs);

double **benchmark(morpho_func_t morphos[], p_struct_elem_dim s,long nb_funcs, const long nb_tests, long min_size, long max_size, long step, const char *filename_prefix, int save_output);

unsigned long long get_cpu_cycles(morpho_func_t morpho, uint8 **ppInput,  long nrl, long nrh, long ncl, long nch, p_struct_elem_dim s, uint8 **ppOutput);
unsigned long long get_min_cycle(unsigned long long *cycles, long packet_size);

uint8 Morpho_Test_5x5_Rect(morpho_func_t morpho, uint8 **ppInput, p_struct_elem_dim s);
void test_morpho(morpho_func_t erosion, morpho_func_t dilation);
#endif /*  __TEST_MORPHO_H__ */