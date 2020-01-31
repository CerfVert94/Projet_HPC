#ifndef  __TEST_BENCHMARK_H__
#define __TEST_BENCHMARK_H__

#pragma message("  include  benchmark.h")

unsigned long long get_cpu_cycles_of_morpho           (struct morpho_set *ptr_mset                  , uint8 **X,  long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y);
unsigned long long get_min_cpu_cycles_of_morpho       (struct morpho_set *ptr_mset, long packet_size, uint8 **X,  long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y);

unsigned long long get_cpu_cycles_of_packed_morpho    (struct morpho_set *ptr_mset                  , uint8 **X,  long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y);
unsigned long long get_min_cpu_cycles_of_packed_morpho(struct morpho_set *ptr_mset, long packet_size, uint8 **X,  long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y);

unsigned long long get_cpu_cycles_of_sd               (struct morpho_set *ptr_mset                  , uint8** X, uint8** Y, uint8** Z, long nrl, long nrh, long ncl, long nch);
unsigned long long get_min_cpu_cycles_sd             (struct morpho_set *ptr_mset, long packet_size, uint8** X, uint8** Y, uint8** Z, long nrl, long nrh, long ncl, long nch);

double **benchmark_of_morpho        (struct morpho_set *msets, long nb_sets, long ls, long hs, long step, int nb_tests, int packet_size);
double **benchmark_of_packed_morpho (struct morpho_set *msets, long nb_sets, long ls, long hs, long step, int nb_tests, int packet_size);
double **benchmark_of_sd            (struct morpho_set *msets, long nb_sets, long ls, long hs, long step, int nb_tests, int packet_size);

double **init_benchmark_results(long nb_funcs, long size, long step);
void     free_benchmark_results(double **results, long nb_funcs);






#endif /*_BENCHMARK_H_*/