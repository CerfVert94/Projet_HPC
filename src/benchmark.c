
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "nrdef.h"
#include "vnrdef.h"
#include "nrutil.h"
#include "vnrutil.h"
#include "mynrutil.h"
#include "myvnrutil.h"
#include "util.h"
#include "img.h"
#include "img_SIMD.h"
#include <morpho.h>
#include <mouvement.h>
#include <test_morpho.h>
#include <benchmark.h>
#include <mutil.h>
#include <string.h>
#include <assert.h>
#include <x86intrin.h>
#include <time.h>
#include <math.h>


void create_false_inputs(p_image *t0, p_image *t1, uint8 ***X, uint8 ***Y, uint8 ***Z, p_vimage *v_t0, p_vimage *v_t1, vuint8 ***vX, vuint8 ***vY, vuint8 ***vZ, long size, int *i0, int *i1, int *v0, int *v1);
void free_false_inputs(p_image t0, p_image t1, uint8 **X, uint8 **Y, uint8 **Z, p_vimage v_t0, p_vimage v_t1, vuint8 **vX, vuint8 **vY, vuint8 **vZ, long size, int i0, int i1, int v0, int v1);

void create_false_inputs(p_image *t0, p_image *t1, uint8 ***X, uint8 ***Y, uint8 ***Z, p_vimage *v_t0, p_vimage *v_t1, vuint8 ***vX, vuint8 ***vY, vuint8 ***vZ, long size, int *i0, int *i1, int *v0, int *v1)
{
	
	*X    = ui8matrix_checker(0 - BORD, size + BORD, 0 - BORD, size + BORD, 3, 1); 
	*Y    = ui8matrix_checker(0 - BORD, size + BORD, 0 - BORD, size + BORD, 3, 0); 
	*Z    = ui8matrix        (0 - BORD, size + BORD, 0 - BORD, size + BORD); 
	
	*vX   = ui8matrix_to_vui8matrix(*X, 0 - BORD, size + BORD, 0 - BORD, size + BORD, i0, i1, v0, v1);
	*vY   = ui8matrix_to_vui8matrix(*Y, 0 - BORD, size + BORD, 0 - BORD, size + BORD, i0, i1, v0, v1);
	*vZ   =              vui8matrix(    0 - BORD, size + BORD, *v0, *v1);
	*t0   = create_image_from_ui8matrix( *X, 0, size, 0, size);
	*t1   = create_image_from_ui8matrix( *Y, 0, size, 0, size);
	*v_t0 = create_vimage_from_ui8matrix(*X, 0, size, 0, size);
	*v_t1 = create_vimage_from_ui8matrix(*Y, 0, size, 0, size);
	// printf("create: %d %d %d %d\n", *i0,*i1,*v0,*v1);
	// print_vui8matrix((*v_t0)->I, (*v_t0)->nrl, (*v_t0)->nrh, (*v_t0)->v0, (*v_t0)->v1, "%u", "Image");getchar();
	

}
void free_false_inputs(p_image t0, p_image t1, uint8 **X, uint8 **Y, uint8 **Z, p_vimage v_t0, p_vimage v_t1, vuint8 **vX, vuint8 **vY, vuint8 **vZ, long size, int i0, int i1, int v0, int v1)
{
	free_ui8matrix(  X, 0 - BORD, size + BORD, 0 - BORD, size + BORD);
	free_ui8matrix(  Y, 0 - BORD, size + BORD, 0 - BORD, size + BORD);
	free_ui8matrix(  Z, 0 - BORD, size + BORD, 0 - BORD, size + BORD);
	// printf("create: %d %d %d %d\n", i0,i1,v0,v1);
	free_vui8matrix(vX, i0, i1, v0, v1);
	free_vui8matrix(vY, i0, i1, v0, v1);
	free_vui8matrix(vZ, i0, i1, v0, v1);
	free_image(t0);	
	free_image(t1);
	free_vimage(v_t0);
	free_vimage(v_t1);
}

void save_benchmark(const char *filename, void *sets, size_t struct_size, int nb_sets, double **results, long min_size, long max_size, long step);

void launch_complete_process_benchmark(const char *filename, struct complete_process_set *cps, const int nb_sets, const int nb_tests,const int packet_size, long min_size, long max_size, long step) {
    double **results;
    results = benchmark_of_complete_process(cps, nb_sets, min_size, max_size, step, nb_tests, packet_size);
    save_benchmark(filename, cps, sizeof(struct complete_process_set), nb_sets, results, min_size, max_size, step);
}

void launch_morpho_benchmark(const char *filename, struct morpho_set *morphos, const int nb_sets, const int nb_tests,const int packet_size, long min_size, long max_size, long step) {
    double **results;
    results = benchmark_of_morpho(morphos, nb_sets, min_size, max_size, step, nb_tests, packet_size);
    save_benchmark(filename, morphos, sizeof(struct morpho_set), nb_sets, results, min_size, max_size, step);
}


void launch_SD_step_benchmark(const char *filename, struct sd_set *sd_steps, const int nb_sets, const int nb_tests,const int packet_size, long min_size, long max_size, long step) {
    double **results;
    results = benchmark_of_sd_step(sd_steps, nb_sets, min_size, max_size, step, nb_tests, packet_size);
    save_benchmark(filename, sd_steps, sizeof(struct sd_set), nb_sets, results, min_size, max_size, step);
}


void launch_SD_benchmark(const char *filename, struct complete_sd_set *csds, const int nb_sets, const int nb_tests,const int packet_size, long min_size, long max_size, long step) {
    double **results;
    results = benchmark_of_sd(csds, nb_sets, min_size, max_size, step, nb_tests, packet_size);
	save_benchmark(filename, csds, sizeof(struct complete_sd_set), nb_sets, results, min_size, max_size, step);
}

void launch_packed_morpho_benchmark(const char *filename, struct morpho_set *packed_morphos, const int nb_sets, const int nb_tests,const int packet_size, long min_size, long max_size, long step) {
    // double **results;
    // results = benchmark_o_(morphos, nb_sets, min_size, max_size, step, nb_tests, packet_size);
    // save_benchmark(filename, morphos, sizeof(struct morpho_set), nb_sets, results, min_size, max_size, step);
}
double **init_benchmark_results(long nb_funcs, long size, long step)
{
	double **results;
	results = (double **) malloc(sizeof(double *) * nb_funcs);
	if (!results)    
		exit_on_error("malloc failed");

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
unsigned long long get_min_cycles(unsigned long long cycles[],long packet_size)
{
	unsigned long long min_cycles = cycles[0];

	// printf("cycles0 : %llu\n",cycles[0]);
	for (long i = 1; i < packet_size; i++) {

	// printf("cycles%d : %llu\n",i,cycles[i]);
		if (min_cycles > cycles[i]) 
			min_cycles = cycles[i];
	}
	// printf("min : %llu\n",min_cycles);
	return min_cycles;
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
	min_cycles = get_min_cycles(cycles, packet_size);
	free(cycles);
	return min_cycles;
}
unsigned long long get_cpu_cycles_of_sd(struct complete_sd_set *csdset, p_image t0, p_image t1, uint8 n_coeff, uint8 v_min, uint8 v_max)
{
	unsigned long long begin = 0, end = 0, cycles = 0;

	begin = __rdtsc();
	csdset->sd_step0(t0->M, t0->I, t0->V, t0->nrl, t0->nrh, t0->ncl, t0->nch, n_coeff, v_min, v_max);
	end = __rdtsc() - begin;
	begin = __rdtsc() + end ;
	csdset->sd_func(t0, t1, n_coeff, v_min, v_max);
	end = __rdtsc();
	return end - begin;
}
unsigned long long get_min_cpu_cycles_of_sd(struct complete_sd_set *csdset, long packet_size, p_image t0, p_image t1, uint8 n_coeff, uint8 v_min, uint8 v_max)
{
	unsigned long long min_cycles, *cycles;
	int i;
	if(packet_size <= 0) {
		fprintf(stderr, "Error : the size of a packet should be at least 1.\n");
		exit(EXIT_FAILURE);
	}
	cycles = (unsigned long long *) malloc(sizeof(unsigned long long) * packet_size);
	for (int i = 0; i < packet_size; i++) 
		cycles[i] = get_cpu_cycles_of_sd(csdset, t0, t1, n_coeff, v_min, v_max);

	min_cycles = get_min_cycles(cycles, packet_size);
	free(cycles);
	return min_cycles;
}

unsigned long long get_cpu_cycles_of_sd_step(struct sd_set *sdset, uint8** X, uint8** Y, uint8** Z, long nrl, long nrh, long ncl, long nch, uint8 n_coeff, uint8 v_min, uint8 v_max)
{
    unsigned long long begin = 0, end = 0, cycles = 0;

	begin = __rdtsc();
	sdset->sd_func(X, Y, Z, nrl, nrh, ncl, nch, n_coeff, v_min, v_max);
	end = __rdtsc();
	return end - begin;
}

unsigned long long get_min_cpu_cycles_of_sd_step(struct sd_set *sdset, long packet_size, uint8** X, uint8** Y, uint8** Z, long nrl, long nrh, long ncl, long nch, uint8 n_coeff, uint8 v_min, uint8 v_max)
{
	unsigned long long min_cycles, *cycles;
	int i;
	if(packet_size <= 0) {
		fprintf(stderr, "Error : the size of a packet should be at least 1.\n");
		exit(EXIT_FAILURE);
	}
	cycles = (unsigned long long *) malloc(sizeof(unsigned long long) * packet_size);
	for (int i = 0; i < packet_size; i++) 
		cycles[i] = get_cpu_cycles_of_sd_step(sdset, X, Y, Z, nrl, nrh, ncl, nch, n_coeff, v_min, v_max);

	min_cycles = get_min_cycles(cycles, packet_size);
	free(cycles);
	return min_cycles;
}



unsigned long long get_cpu_cycles_of_vec_morpho(struct morpho_set *ptr_mset, vuint8 **vX, int i0, int i1, long ncl, long nch, int j0, int j1, vuint8 **vTempBuffer, vuint8 **vY)
{
    unsigned long long begin = 0, end = 0, cycles = 0;
	begin = __rdtsc();
	ptr_mset->vec_morpho_func(vX, i0, i1, ncl, nch, j0, j1, vTempBuffer, vY);
	end = __rdtsc();
	return end - begin;
}

unsigned long long get_min_cpu_cycles_of_vec_morpho(struct morpho_set *ptr_mset, long packet_size, vuint8 **vX, int i0, int i1, long ncl, long nch, int j0, int j1, vuint8 **vTempBuffer, vuint8 **vY)
{
	unsigned long long min_cycles, *cycles;
	int i;
	if(packet_size <= 0) {
		fprintf(stderr, "Error : the size of a packet should be at least 1.\n");
		exit(EXIT_FAILURE);
	}
	cycles = (unsigned long long *) malloc(sizeof(unsigned long long) * packet_size);
	for (int i = 0; i < packet_size; i++) {
		cycles[i] = get_cpu_cycles_of_vec_morpho(ptr_mset, vX, i0, i1, ncl, nch, j0, j1, vTempBuffer, vY);
	}
	min_cycles = get_min_cycles(cycles, packet_size);
	free(cycles);
	return min_cycles;
}
unsigned long long get_cpu_cycles_of_vec_sd(struct complete_sd_set *csdset, p_vimage t0, p_vimage t1, uint8 n_coeff, uint8 v_min, uint8 v_max)
{
	unsigned long long begin = 0, end = 0, cycles = 0;

	begin = __rdtsc();
	csdset->vec_sd_step0(t0->M, t0->I, t0->V, t0->nrh, t0->nrh, t0->v0, t0->v1, n_coeff, v_min, v_max);
	end = __rdtsc() - begin;
	begin = __rdtsc() + end ;
	csdset->vec_sd_func(t0, t1, n_coeff, v_min, v_max);
	end = __rdtsc();
	return end - begin;
}


unsigned long long get_min_cpu_cycles_of_vec_sd(struct complete_sd_set *csdset, long packet_size, p_vimage t0, p_vimage t1, uint8 n_coeff, uint8 v_min, uint8 v_max)
{
	unsigned long long min_cycles, *cycles;
	int i;
	if(packet_size <= 0) {
		fprintf(stderr, "Error : the size of a packet should be at least 1.\n");
		exit(EXIT_FAILURE);
	}
	cycles = (unsigned long long *) malloc(sizeof(unsigned long long) * packet_size);
	for (int i = 0; i < packet_size; i++){
		cycles[i] = get_cpu_cycles_of_vec_sd(csdset, t0, t1, n_coeff, v_min, v_max);
		
	}

	min_cycles = get_min_cycles(cycles, packet_size);
	free(cycles);
	return min_cycles;
}

unsigned long long get_cpu_cycles_of_vec_sd_step(struct sd_set *sdset, vuint8 **vX, vuint8 **vY, vuint8 **vZ, int i0, int i1, int j0, int j1, uint8 n_coeff, uint8 v_min, uint8 v_max)
{
    unsigned long long begin = 0, end = 0, cycles = 0;

	begin = __rdtsc();
	sdset->vec_sd_func(vX, vY, vZ, i0, i1, j0, j1, n_coeff, v_min, v_max);
	end = __rdtsc();
	return end - begin;
}

unsigned long long get_min_cpu_cycles_of_vec_sd_step(struct sd_set *sdset, long packet_size, vuint8 **vX, vuint8 **vY, vuint8 **vZ, int i0, int i1, int j0, int j1, uint8 n_coeff, uint8 v_min, uint8 v_max)
{
	unsigned long long min_cycles, *cycles;
	int i;
	if(packet_size <= 0) {
		fprintf(stderr, "Error : the size of a packet should be at least 1.\n");
		exit(EXIT_FAILURE);
	}
	cycles = (unsigned long long *) malloc(sizeof(unsigned long long) * packet_size);
	for (int i = 0; i < packet_size; i++) 
		cycles[i] = get_cpu_cycles_of_vec_sd_step(sdset, vX, vY, vZ, i0, i1, j0, j1, n_coeff, v_min, v_max);

	min_cycles = get_min_cycles(cycles, packet_size);
	free(cycles);
	return min_cycles;
}

unsigned long long get_cpu_cycles_of_complete_process(struct complete_process_set *cproc)
{
    unsigned long long begin = 0, end = 0, cycles = 0;	
	uint8 **tempBuffer;
    p_image t0, t1;

    tempBuffer = cproc->Y;
    t0 = cproc->t0;
	t1 = cproc->t1;

	
	begin = __rdtsc() + end;
	cproc->sd_step0(t0->M, t0->I, t0->V, t0->nrl, t0->nrh, t0->ncl, t0->nch, N, Vmin, Vmax);
    cproc->sd_func(t0, t1, N, Vmin, Vmax);
    cproc->morpho_func(t1->E, t1->nrl + BORD, t1->nrh - BORD, t1->ncl + BORD, t1->nch - BORD, tempBuffer, t1->Omega);
	end = __rdtsc();
	return end - begin;
}
unsigned long long get_cpu_cycles_of_vec_complete_process(struct complete_process_set *cproc)
{
    unsigned long long begin = 0, end = 0, cycles = 0;	
	vuint8 **vTempBuffer;
    p_vimage t0, t1;

    vTempBuffer = cproc->vY;
    t0 = cproc->v_t0;
	t1 = cproc->v_t1;
	
	begin = __rdtsc() + end;
    cproc->vec_sd_step0(t0->M, t0->I, t0->V, t0->nrl, t0->nrh, t0->v0, t0->v1, N, Vmin, Vmax);
	cproc->vec_sd_func(t0, t1, N, Vmin, Vmax);
    cproc->vec_morpho_func(t1->E, (int)t1->nrl + BORD, (int)t1->nrh - BORD, (int)t1->ncl + BORD, (int)t1->nch - BORD, t1->v0 + vBORD, t1->v1 - vBORD, vTempBuffer, t1->Omega);
	end = __rdtsc();
	return end - begin;
}


unsigned long long get_min_cpu_cycles_of_vec_complete_process(struct complete_process_set *cproc, long packet_size)
{
	unsigned long long min_cycles, *cycles;
	int i;

	if(packet_size <= 0) {
		fprintf(stderr, "Error : the size of a packet should be at least 1.\n");
		exit(EXIT_FAILURE);
	}
	
	cycles = (unsigned long long *) malloc(sizeof(unsigned long long) * packet_size);

	
	for (int i = 0; i < packet_size; i++) 
		cycles[i] = get_cpu_cycles_of_vec_complete_process(cproc);

	min_cycles = get_min_cycles(cycles, packet_size);
	free(cycles);
	return min_cycles;
}
unsigned long long get_min_cpu_cycles_of_complete_process(struct complete_process_set *cproc, long packet_size)
{
	unsigned long long min_cycles, *cycles;
	int i;

	if(packet_size <= 0) {
		fprintf(stderr, "Error : the size of a packet should be at least 1.\n");
		exit(EXIT_FAILURE);
	}

	cycles = (unsigned long long *) malloc(sizeof(unsigned long long) * packet_size);

	
	for (int i = 0; i < packet_size; i++) 
		cycles[i] = get_cpu_cycles_of_complete_process(cproc);

	min_cycles = get_min_cycles(cycles, packet_size);
	free(cycles);
	return min_cycles;
}


double **benchmark_of_complete_process(struct complete_process_set *cproc, long nb_sets, long ls, long hs, long step, int nb_tests, int packet_size)
{
	unsigned long long  min_cycles_sum, begin, end;
	long size, i, idx_test, nrl, nrh, ncl, nch, cnt = 0;
	double **results;
	// Initialize to save the benchmark result.
	results = init_benchmark_results(nb_sets, (hs - ls)  + 1, step);
	
	cnt = 0;
	for (size = ls - 1; size < hs; size += step) {
		// idx_set = 0;
		for (i = 0; i < nb_sets; i++) 
		{
			create_false_inputs(&cproc[i].t0  , &cproc[i].t1  , &cproc[i].X , &cproc[i].Y , &cproc[i].Z, 
			                    &cproc[i].v_t0, &cproc[i].v_t1, &cproc[i].vX, &cproc[i].vY, &cproc[i].vZ, 
								size, &cproc[i].i0, &cproc[i].i1, &cproc[i].j0, &cproc[i].j1);

			
			
			cproc[i].nrl = 0; cproc[i].nrh = size;
			cproc[i].ncl = 0; cproc[i].nch = size;
	
			if (cproc[i].instr_type == SCALAR) {
				// p_image t0 = cproc[i].t0;
				begin = __rdtsc();			
				min_cycles_sum = 0;
				for (idx_test = 0; idx_test < nb_tests; idx_test++)
					min_cycles_sum += get_min_cpu_cycles_of_complete_process(&cproc[i], packet_size);
				end = __rdtsc();
			}
			else if (cproc[i].instr_type == SIMD) {
				// p_vimage t0 = cproc[i].v_t0;
				begin = __rdtsc();			
				min_cycles_sum = 0;
				for (idx_test = 0; idx_test < nb_tests; idx_test++)
					min_cycles_sum += get_min_cpu_cycles_of_vec_complete_process(&cproc[i], packet_size);
				end = __rdtsc();
			}
			// printf("%llu / %llu\n", min_cycles_sum, (nb_tests * (size + 1) * (size + 1)));
			results[i][cnt] = ((double)min_cycles_sum / (nb_tests * (size + 1) * (size + 1)));
			
			
			
			if ((size + 1) % 500 == 0 || size >= hs - 1) 
				printf("\t["LALIGNED_STR"] Ran SigmaDelta %d * %d times on %ld x %ld matrix during %llu cycles.\n",  cproc[i].func_name, packet_size, nb_tests, size + 1, size + 1, (end - begin));			

			free_false_inputs(cproc[i].t0  , cproc[i].t1  , cproc[i].X , cproc[i].Y , cproc[i].Z, 
			                  cproc[i].v_t0, cproc[i].v_t1, cproc[i].vX, cproc[i].vY, cproc[i].vZ, 
							  size, cproc[i].i0, cproc[i].i1, cproc[i].j0, cproc[i].j1);
		}
		cnt++;
	}
	return results;
}

double **benchmark_of_sd_step(struct sd_set     *sdsets, long nb_sets, long ls, long hs, long step, int nb_tests, int packet_size)
{
	unsigned long long  min_cycles_sum, begin, end;
	long size, idx_set, idx_test, nrl, nrh, ncl, nch, cnt = 0;
	double **results;
	uint8 **X, **Y, **Z;
	vuint8 **vX, **vY, **vZ;
	uint8 n_coeff = 2;
	uint8 v_min = 1;
	uint8 v_max = 254;
	int i0, i1, v0, v1;
	// Initialize to save the benchmark result.
	results = init_benchmark_results(nb_sets, (hs - ls)  + 1, step);
	
	cnt = 0;
	for (size = ls - 1; size < hs; size += step) {
		// idx_set = 0;
		for (idx_set = 0; idx_set < nb_sets; idx_set++) 
		{
			X = ui8matrix_checker(0, size, 0, size, 3, 1); 
			Y = ui8matrix_checker(0, size, 0, size, 3, 0); 
			Z = ui8matrix        (0, size, 0, size); 
			vX = ui8matrix_to_vui8matrix(X, 0, size, 0, size, &i0, &i1, &v0, &v1);
			vY = ui8matrix_to_vui8matrix(Y, 0, size, 0, size, &i0, &i1, &v0, &v1);
			vZ =              vui8matrix(   0, size, 0, size);
			n_coeff = sdsets[idx_set].n_coeff;
			v_min   = sdsets[idx_set].v_min;
			v_max   = sdsets[idx_set].v_max;
			
			if (sdsets[idx_set].instr_type == SCALAR) {
				begin = __rdtsc();			
				min_cycles_sum = 0;
				for (idx_test = 0; idx_test < nb_tests; idx_test++)
					min_cycles_sum += get_min_cpu_cycles_of_sd_step(&sdsets[idx_set], packet_size, X, Y, Z, 0, size, 0, size, n_coeff, v_min, v_max);
				results[idx_set][cnt] = ((double)min_cycles_sum / (nb_tests * (size + 1) * (size + 1)));
				end = __rdtsc();
			}
			else if (sdsets[idx_set].instr_type == SIMD) {
				begin = __rdtsc();			
				min_cycles_sum = 0;
				for (idx_test = 0; idx_test < nb_tests; idx_test++)
					min_cycles_sum += get_min_cpu_cycles_of_vec_sd_step(&sdsets[idx_set], packet_size, vX, vY, vZ, i0, i1, v0, v1, n_coeff, v_min, v_max);
				results[idx_set][cnt] = ((double)min_cycles_sum / (nb_tests * (size + 1) * (size + 1)));
				end = __rdtsc();
			}
			
			if ((size + 1) % 500 == 0 || size >= hs - 1) 
				printf("\t["LALIGNED_STR"] Ran SigmaDelta %d * %d times on %ld x %ld matrix during %llu cycles.\n",  sdsets[idx_set].func_name, packet_size, nb_tests, size + 1, size + 1, (end - begin));			

			free_ui8matrix(  X, 0, size, 0, size);
			free_ui8matrix(  Y, 0, size, 0, size);
			free_ui8matrix(  Z, 0, size, 0, size);
			free_vui8matrix(vX, 0, size, 0, size);
			free_vui8matrix(vY, 0, size, 0, size);
			free_vui8matrix(vZ, 0, size, 0, size);
		}
		
		cnt++;
	}
	return results;
}

double **benchmark_of_morpho(struct morpho_set *morphos, long nb_sets, long ls, long hs, long step, int nb_tests, int packet_size)
{
	unsigned long long  min_cycles_sum, begin, end;
	long size, idx_set, idx_test, nrl, nrh, ncl, nch, cnt = 0;
	double **results;
	uint8 **X, **Y, **temp_buffer;
	vuint8 **vX, **vY, **vTempBuffer;
	int i0, i1, j0, j1;
	int k0, k1, l0, l1;
	
	// Initialize to save the benchmark result.
	results = init_benchmark_results(nb_sets, (hs - ls)  + 1, step);
	
	cnt = 0;
	for (size = ls - 1; size < hs; size += step) {
		
		for (idx_set = 0; idx_set < nb_sets; idx_set++) {
			X            = ui8matrix_checker(-2, size + 2, -2, size + 2, 3, 1); 
			Y            = ui8matrix(		 -2, size + 2, -2, size + 2);
			temp_buffer  = ui8matrix(		 -2, size + 2, -2, size + 2);
			vX			 = ui8matrix_to_vui8matrix(X, -2, size + 2, -2, size + 2, &i0, &i1, &j0, &j1);
			vY			 = ui8matrix_to_vui8matrix(Y, -2, size + 2, -2, size + 2, &k0, &k1, &l0, &l1);
			vTempBuffer  = vui8matrix(i0, i1, j0, j1); 
			if (morphos[idx_set].instr_type == SCALAR) {
				begin = __rdtsc();			
				min_cycles_sum = 0;
				for (idx_test = 0; idx_test < nb_tests; idx_test++)
					min_cycles_sum += get_min_cpu_cycles_of_morpho(&morphos[idx_set], packet_size, X, 0, size, 0, size, temp_buffer, Y);
				end = __rdtsc();			
			}
			else if (morphos[idx_set].instr_type == SIMD) {
				begin = __rdtsc();			
				min_cycles_sum = 0;
				for (idx_test = 0; idx_test < nb_tests; idx_test++)
					min_cycles_sum += get_min_cpu_cycles_of_vec_morpho(&morphos[idx_set], packet_size, vX, 0, size, 0, size, j0 + 1, j1, vTempBuffer, vY);
				end = __rdtsc();	
			}

			results[idx_set][cnt] = ((double)min_cycles_sum / (nb_tests * (size + 1) * (size + 1)));
			if ((size + 1) % 500 == 0 || size >= hs - 1) 
				printf("\t["LALIGNED_STR"] Ran morpho %d * %d times on %ld x %ld matrix during %llu cycles (min : %2.02lf).\n",  morphos[idx_set].func_name, packet_size, nb_tests, size + 1, size + 1, (end - begin), results[idx_set][cnt]);			
			free_vui8matrix(vTempBuffer, i0, i1, j0, j1);
			free_vui8matrix(vX, i0, i1, j0, j1);
			free_vui8matrix(vY, k0, k1, l0, l1);
			free_ui8matrix(temp_buffer, -2, size + 2, -2, size + 2);
			free_ui8matrix(X, -2, size + 2, -2, size + 2);
			free_ui8matrix(Y, -2, size + 2, -2, size + 2);
		} 

		cnt++;
	}
	return results;
}


double **benchmark_of_packed_morpho(struct morpho_set *morphos, long nb_sets, long ls, long hs, long step, int nb_tests, int packet_size)
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
				min_cycles_sum += get_min_cpu_cycles_of_morpho(&morphos[idx_set], packet_size, packedX, packed_nrl, packed_nrh, packed_ncl, packed_nch, temp_buffer, Y);
			
			results[idx_set][cnt] = ((double)min_cycles_sum / (nb_tests * (size + 1) * (size + 1)));
			unfcpack_ui8matrix_ui8matrix(Y, 0, size, 0, size, packed_nrl, packed_nrh, packed_ncl, packed_nch, bord, Z);
			end = __rdtsc();
			if ((size + 1) % 500 == 0 || size >= hs - 1) 
				printf("\t["LALIGNED_STR"] Ran morpho %d * %d times on %ld x %ld matrix during %llu cycles (min : %2.02lf).\n",  morphos[idx_set].func_name, packet_size, nb_tests, size + 1, size + 1, (end - begin), results[idx_set][cnt]);			
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


double **benchmark_of_sd(struct complete_sd_set *csdsets, long nb_sets, long ls, long hs, long step, int nb_tests, int packet_size)
{
	unsigned long long  min_cycles_sum, begin, end;
	long size, idx_set, idx_test, nrl, nrh, ncl, nch, cnt = 0;
	double **results;
	uint8 **X, **Y;
	// vuint8 **vX, **vY;
	p_image t0, t1;
	p_vimage vec_t0, vec_t1;
	
	int n_coeff = 2;
	int v_min = 1;
	int v_max = 254;
	int i0, i1, v0, v1;
	
	// Initialize to save the benchmark result.
	results = init_benchmark_results(nb_sets, (hs - ls)  + 1, step);
	
	cnt = 0;
	for (size = ls - 1; size < hs; size += step) {
		for (idx_set = 0; idx_set < nb_sets; idx_set++) {

			X = ui8matrix_checker(0, size, 0, size, 3, 1); 
			Y = ui8matrix_checker(0, size, 0, size, 3, 0); 
			// vX = ui8matrix_to_vui8matrix(X, 0, size, 0, size, &i0, &i1, &v0, &v1);
			// vY = ui8matrix_to_vui8matrix(Y, 0, size, 0, size, &i0, &i1, &v0, &v1);
			
			t0 = create_image_from_ui8matrix(X, 0, size, 0, size);
			t1 = create_image_from_ui8matrix(Y, 0, size, 0, size);

			n_coeff = csdsets[idx_set].n_coeff;
			v_min   = csdsets[idx_set].v_min;
			v_max   = csdsets[idx_set].v_max;
			
			if (csdsets[idx_set].instr_type == SCALAR) {
				begin = __rdtsc();			
				min_cycles_sum = 0;
				for (idx_test = 0; idx_test < nb_tests; idx_test++)
					min_cycles_sum += get_min_cpu_cycles_of_sd(&csdsets[idx_set], packet_size, t0, t1, n_coeff, v_min, v_max);
				end = __rdtsc();
			}
			else if (csdsets[idx_set].instr_type == SIMD) {
				vec_t0 = create_vimage_from_ui8matrix(X, 0, size, 0, size);
				vec_t1 = create_vimage_from_ui8matrix(Y, 0, size, 0, size);
				begin = __rdtsc();			
				min_cycles_sum = 0;
				for (idx_test = 0; idx_test < nb_tests; idx_test++);
					min_cycles_sum += get_min_cpu_cycles_of_vec_sd(&csdsets[idx_set], packet_size, vec_t0, vec_t1, n_coeff, v_min, v_max);
				end = __rdtsc();
				free_vimage(vec_t0);
				free_vimage(vec_t1);
			}
			results[idx_set][cnt] = ((double)min_cycles_sum / (nb_tests * (size + 1) * (size + 1)));
			if ((size + 1) % 500 == 0 || size >= hs - 1) 
				printf("\t["LALIGNED_STR"] Ran SigmaDelta %d * %d times on %ld x %ld matrix during %llu cycles (min : %2.02lf).\n",  csdsets[idx_set].func_name, packet_size, nb_tests, size + 1, size + 1, (end - begin), results[idx_set][cnt]);			

			free_image(t0);
			free_image(t1);
			free_ui8matrix(X, 0, size, 0, size);
			free_ui8matrix(Y, 0, size, 0, size);
		}
		cnt++;
	}
	return results;
}
void save_benchmark(const char *filename, void *sets, size_t struct_size, int nb_sets, double **results, long min_size, long max_size, long step)
{
	FILE *file = fopen(filename, "w");
    if (!file)
        exit_on_error("fopen failed");
    long k = 0;
	char* set_ptr = (char*)sets;

    fprintf(file, "#%*s ", 4, "Size");
    for (long i = 0; i < nb_sets; i++) {
		
		// printf(RALIGNED_STR, (char*)&set_ptr[i * struct_size]);
        fprintf(file, RALIGNED_STR" ", (char*)&set_ptr[i * struct_size]);
    }

    fputc('\n', file);
    for (long size = min_size; size < max_size + 1; size += step) {
        fprintf(file, "%*ld ", 4, size);
        for (long i = 0; i < nb_sets; i++) {
            fprintf(file, RALIGNED_DOUBLE" ", results[i][k]);
         }
        fprintf(file, "\n");
        k++;
    }
    free_benchmark_results(results, 1);
    fclose(file);
}
