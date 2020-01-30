#include <nrdef.h>
#include <nrutil.h>
#include <img.h>
#include <morpho.h>
#include <stdbool.h>
#include <test_morpho.h>
#include <benchmark.h>
#include <stdio.h>
#include <stdlib.h>
#include <util.h>
#include <mutil.h>
#include <string.h>
#include <assert.h>
#include <x86intrin.h>
#include <time.h>
#include <math.h>


double **init_benchmark_results(long nb_funcs, long size, long step)
{
	double **results;
	results = (double **) malloc(sizeof(double *) * nb_funcs);
	if (!results)    
		exit_on_error("malloc failed");

	printf("Nb : %ld\n", (size / step) + 1); 
	results[0] = (double *) malloc(sizeof(double) * nb_funcs * (size / step) + 1);
	if (!results[0]) 
		exit_on_error("malloc failed");
	for (long i = 1; i < nb_funcs; i++) results[i] = results[i - 1] + (size / step);
	return results;

}
void free_benchmark_results(double **results, long nb_funcs)
{
	free((char *)results[0]);
	free((char **)results);
}

unsigned long long get_cpu_cycles_of_morpho(struct morpho_set *ptr_mset, uint8 **X,  long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
    unsigned long long begin = 0, end = 0, cycles = 0;

	begin = __rdtsc();
	ptr_mset->morpho_func(X, nrl, nrh, ncl, nch, temp_buffer, Y);
	end = __rdtsc();
	return end - begin;
}

unsigned long long get_min_cpu_cycles_of_morpho(struct morpho_set *ptr_mset, long packet_size, uint8 **X,  long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y)
{
	unsigned long long min_cycles, *cycles;
	int i;
	if(packet_size <= 0) {
		fprintf(stderr, "Error : the size of a packet should be at least 1.\n");
		exit(EXIT_FAILURE);
	}
	cycles = (unsigned long long *) malloc(sizeof(unsigned long long) * packet_size);
	for (int i = 0; i < packet_size; i++) {
		cycles[i] = get_cpu_cycles_of_morpho(ptr_mset, X, nrl, nrh, ncl, nch, temp_buffer, Y);
	}

	min_cycles = cycles[0];
	for (long i = 1; i < packet_size; i++) 
		if (min_cycles > cycles[i]) 
			min_cycles = cycles[i];
	free(cycles);
	return min_cycles;
}





double **benchmark_of_morpho(struct morpho_set *msets, long nb_sets, long ls, long hs, long step, int nb_tests, int packet_size)
{
	unsigned long long  min_cycles_sum, begin, end;
	long size, idx_set, idx_test, nrl, nrh, ncl, nch, cnt = 0;
	double **results;
	uint8 **X, **Y, **temp_buffer;
	
	// Initialize to save the benchmark result.
	results = init_benchmark_results(nb_sets, (hs - ls)  + 1, step);
	// the least largest 3x3 checkered square matrix as an input / output matrix.
	
	cnt = 0;
	for (size = ls - 1; size < hs; size += step) {
		X            = ui8matrix_checker(-2, size + 2, -2, size + 2, 3, 1); 
		Y            = ui8matrix(0, size, 0, size);
		temp_buffer  = ui8matrix(-2, size + 2, -2, size + 2); 
		
		for (idx_set = 0; idx_set < nb_sets; idx_set++) {
			begin = __rdtsc();			
			// memset_ui8matrix(Y, 0, 0, size, 0, size);

	
			// printf("size:%ld / idx : %d\r\n",size,idx_set);
			// if(idx_set == 7){
			// 	getchar();
			// }
			min_cycles_sum = 0;
			for (idx_test = 0; idx_test < nb_tests; idx_test++)
				min_cycles_sum += get_min_cpu_cycles_of_morpho(&msets[idx_set], packet_size, X, 0, size, 0, size, temp_buffer, Y);
			
			results[idx_set][cnt] = ((double)min_cycles_sum / (nb_tests * (size + 1) * (size + 1)));
			end = __rdtsc();
			if ((size + 1) % 500 == 0 || size >= hs - 1) 
				printf("\t["LALIGNED_STR"] Ran morpho %d * %d times on %ld x %ld matrix during %llu cycles.\n",  msets[idx_set].func_name, packet_size, nb_tests, size + 1, size + 1, (end - begin));			
		}
	
		free_ui8matrix(temp_buffer, -2, size + 2, -2, size + 2);
		free_ui8matrix(X, -2, size + 2, -2, size + 2);
		free_ui8matrix(Y, 0, size, 0, size);

		cnt++;
	}
	return results;
}

double **benchmark_of_packed_morpho(struct morpho_set *msets, long nb_sets, long ls, long hs, long step, int nb_tests, int packet_size)
{
	unsigned long long  min_cycles_sum, begin, end;
	long size, idx_set, idx_test, nrl, nrh, ncl, nch, cnt = 0;
	long packed_nrl, packed_nrh, packed_ncl, packed_nch, bord;
	double **results;
	uint8 **X, **packedX, **Y, **temp_buffer, **Z;
	
	// Initialize to save the benchmark result.
	results = init_benchmark_results(nb_sets, (hs - ls)  + 1, step);
	// the least largest 3x3 checkered square matrix as an input / output matrix.
	
	cnt = 0;
	for (size = ls - 1; size < hs; size += step) {
		packedX      = fcpacked_ui8matrix(0, size, 0, size, &packed_nrl, &packed_nrh, &packed_ncl, &packed_nch, &bord); 
		temp_buffer  = fcpacked_ui8matrix(0, size, 0, size, &packed_nrl, &packed_nrh, &packed_ncl, &packed_nch, &bord); 
		Y            = fcpacked_ui8matrix(0, size, 0, size, &packed_nrl, &packed_nrh, &packed_ncl, &packed_nch, &bord); 
		X            = ui8matrix_checker(-bord, size + bord, -bord, size + bord, 3, 1); 
		Z            = ui8matrix(0, size, 0, size); 
		
		
		for (idx_set = 0; idx_set < nb_sets; idx_set++) {
			memset_ui8matrix(Y, 0, packed_nrl, packed_nrh, packed_ncl, packed_nch);

			begin = __rdtsc();			
			fcpack_ui8matrix_ui8matrix(X, 0, size, 0, size, packed_nrl, packed_nrh, packed_ncl, packed_nch, bord, packedX);
			min_cycles_sum = 0;
			for (idx_test = 0; idx_test < nb_tests; idx_test++)
				min_cycles_sum += get_min_cpu_cycles_of_morpho(&msets[idx_set], packet_size, packedX, packed_nrl, packed_nrh, packed_ncl, packed_nch, temp_buffer, Y);
			
			results[idx_set][cnt] = ((double)min_cycles_sum / (nb_tests * (size + 1) * (size + 1)));
			unfcpack_ui8matrix_ui8matrix(Y, 0, size, 0, size, packed_nrl, packed_nrh, packed_ncl, packed_nch, bord, Z);
			end = __rdtsc();
			if ((size + 1) % 500 == 0 || size >= hs - 1) 
				printf("\t["LALIGNED_STR"] Ran morpho %d * %d times on %ld x %ld matrix during %llu cycles.\n",  msets[idx_set].func_name, packet_size, nb_tests, size + 1, size + 1, (end - begin));			
		}
	
		// free_ui8matrix(temp_buffer, -2, size + 2, -2, size + 2);
		free_packed_ui8matrix(packedX    , packed_nrl, packed_nrh, packed_ncl, packed_nch, bord);
		free_packed_ui8matrix(temp_buffer, packed_nrl, packed_nrh, packed_ncl, packed_nch, bord);
		free_packed_ui8matrix(Y          , packed_nrl, packed_nrh, packed_ncl, packed_nch, bord);
		free_ui8matrix(X, -2, size + 2, -2, size + 2);
		free_ui8matrix(Z, 0, size, 0, size);

		cnt++;
	}
	return results;
}


