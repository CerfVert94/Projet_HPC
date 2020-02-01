#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "nrdef.h"
#include "nrutil.h"
#include "mynrutil.h"
#include "util.h"
#include "img.h"
#include "mouvement.h"
#include <stdbool.h>
#include "test_mouvement.h"
#include <assert.h>

const char *nom_func;


/*---------------------------------------------------*/
void verify_case_SigmaDelta_step1(struct sd_set *sd, int num_case, const char *str_case, uint8 **X, uint8 **Y, uint8 **Z, uint8 x, uint8 y, bool logging)
/*---------------------------------------------------*/
{
	if(logging)
		printf("Case %d : %s\n", num_case, str_case);
	X[0][0] = x;
	Y[0][0] = y;
	sd->sd_func(X, Y, Z, 0, 0, 0, 0, sd->n_coeff, sd->v_min, sd->v_max);
	if(logging)
		printf("\t[M_t0 = %u] [I_t1 = %u] => [M_t1 = %u]\n", X[0][0], Y[0][0], Z[0][0]);
	assert(SD_step1_produces_valid_output(X[0][0], Y[0][0], Z[0][0], logging) == true);
}


/*---------------------------------------------------*/
void verify_case_SigmaDelta_step2(struct sd_set *sd, int num_case, const char *str_case, uint8 **X, uint8 **Y, uint8 **Z, uint8 x, uint8 y, bool logging)
/*---------------------------------------------------*/
{
	if(logging)
		printf("Case %d : %s\n", num_case, str_case);
	X[0][0] = x;
	Y[0][0] = y;
	sd->sd_func(X, Y, Z, 0, 0, 0, 0, sd->n_coeff, sd->v_min, sd->v_max);
	if(logging)
		printf("\t[M_t1 = %u] [I_t1 = %u] => [O_t1 = %u]\n", X[0][0], Y[0][0], Z[0][0]);
	assert(SD_step2_produces_valid_output(X[0][0], Y[0][0], Z[0][0], logging) == true);
}

/*---------------------------------------------------*/
void verify_case_SigmaDelta_step3(struct sd_set *sd, int num_case, const char *str_case, uint8 **V_t0, uint8 **O_t1, uint8 **V_t1, uint8 v_t0, uint8 o_t1, bool logging)
/*---------------------------------------------------*/
{
	if(logging)
		printf("Case %d : %s\n", num_case, str_case);
	V_t0[0][0] = v_t0;
	O_t1[0][0] = o_t1;
	sd->sd_func(V_t0, O_t1, V_t1, 0, 0, 0, 0, sd->n_coeff, sd->v_min, sd->v_max);
	if(logging)
		printf("\t[V_t0 = %u] [O_t1 = %u] => [V_t1 = %u]\n", V_t0[0][0], O_t1[0][0], V_t1[0][0]);
	assert(SD_step3_produces_valid_output(V_t0[0][0], O_t1[0][0], V_t1[0][0], sd->n_coeff, sd->v_min, sd->v_max, logging) == true);
}

/*---------------------------------------------------*/
void verify_case_SigmaDelta_step4(struct sd_set *sd, int num_case, const char *str_case, uint8 **X, uint8 **Y, uint8 **Z, uint8 x, uint8 y, bool logging)
/*---------------------------------------------------*/
{
	if(logging)
		printf("Case %d : %s\n", num_case, str_case);
	X[0][0] = x;
	Y[0][0] = y;
	sd->sd_func(X, Y, Z, 0, 0, 0, 0, sd->n_coeff, sd->v_min, sd->v_max);
	if(logging)
		printf("\t[O_t1 = %u] [V_t1 = %u] => [E_t1 = %u]\n", X[0][0], Y[0][0], Z[0][0]);
	assert(SD_step4_produces_valid_output(X[0][0], Y[0][0], Z[0][0], logging) == true);
}


/*---------------------------------------------------*/
void test_integration_SigmaDelta_step0(char *filename0, char *filename1, struct sd_set *sd, int nb_sets, bool logging)
/*---------------------------------------------------*/
{
	long nrl0, ncl0, nrh0, nch0;
	long nrl1, ncl1, nrh1, nch1;
	long nrl2, ncl2, nrh2, nch2;
    uint8 **X, **Y, **Z;
	uint8 **packedX, **packedY, **temp_packed_buffer, **unpacked;

    X = LoadPGM_ui8matrix(filename0, &nrl0, &nrh0, &ncl0, &nch0);
	Y = LoadPGM_ui8matrix(filename1, &nrl1, &nrh1, &ncl1, &nch1);
	assert(nrl1 == nrl0 && nrh1 == nrh0 && ncl1 == ncl0 && nch1 == nch0);
	Z = ui8matrix(nrl1, nrh1, ncl1, nch1);
	for (int i = 0; i < nb_sets; i++) {
		printf("Integration test for %s\n", sd[i].func_name);
		sd->sd_func(X, Y, Z, nrl1, nrh1, ncl1, nch1, sd[i].n_coeff, sd[i].v_min, sd[i].v_max);

		for (long row = nrl1; row < nrh1 + 1; row++)
			for (long col = nrl1; col < nrh1 + 1; col++){
				if (logging) {
					printf("[row, col] = [%ld, %ld]\n", row, col);
					printf("[M_t0 = %u] [I_t0 = %u] => [V_t0 = %u]\n", X[row][col], Y[row][col], Z[row][col]);
				}
				assert(SD_step0_produces_valid_output(X[row][col], Y[row][col], Z[row][col], logging));
			}
		printf("Test passed.\n");
	}
	free_ui8matrix(X, nrl0, nrh0, ncl0, nch0);
	free_ui8matrix(Y, nrl1, nrh1, ncl1, nch1);
	free_ui8matrix(Z, nrl1, nrh1, ncl1, nch1);
}

/*---------------------------------------------------*/
void test_integration_SigmaDelta_step1(char *filename0, char *filename1, struct sd_set *sd, int nb_sets, bool logging)
/*---------------------------------------------------*/
{
	long nrl0, ncl0, nrh0, nch0;
	long nrl1, ncl1, nrh1, nch1;
	long nrl2, ncl2, nrh2, nch2;
    uint8 **X, **Y, **Z;
	uint8 **packedX, **packedY, **temp_packed_buffer, **unpacked;

    X = LoadPGM_ui8matrix(filename0, &nrl0, &nrh0, &ncl0, &nch0);
	Y = LoadPGM_ui8matrix(filename1, &nrl1, &nrh1, &ncl1, &nch1);
	assert(nrl1 == nrl0 && nrh1 == nrh0 && ncl1 == ncl0 && nch1 == nch0);
	Z = ui8matrix(nrl1, nrh1, ncl1, nch1);
	for (int i = 0; i < nb_sets; i++) {
		printf("Integration test for %s\n", sd[i].func_name);
		sd[i].sd_func(X, Y, Z, nrl1, nrh1, ncl1, nch1, sd[i].n_coeff, sd[i].v_min, sd[i].v_max);

		for (long row = nrl1; row < nrh1 + 1; row++)
			for (long col = nrl1; col < nrh1 + 1; col++){
				if (logging) {
					printf("[row, col] = [%ld, %ld]\n", row, col);
					printf("[M_t0 = %u] [I_t0 = %u] => [M_t1 = %u]\n", X[row][col], Y[row][col], Z[row][col]);
				}
				assert(SD_step1_produces_valid_output(X[row][col], Y[row][col], Z[row][col], logging));
			}
		printf("Test passed.\n");
	}
	free_ui8matrix(X, nrl0, nrh0, ncl0, nch0);
	free_ui8matrix(Y, nrl1, nrh1, ncl1, nch1);
	free_ui8matrix(Z, nrl1, nrh1, ncl1, nch1);
}

/*---------------------------------------------------*/
void test_integration_SigmaDelta_step2(char *filename0, char *filename1, struct sd_set *sd, int nb_sets, bool logging)
/*---------------------------------------------------*/
{
	long nrl0, ncl0, nrh0, nch0;
	long nrl1, ncl1, nrh1, nch1;
	long nrl2, ncl2, nrh2, nch2;
    uint8 **X, **Y, **Z;
	uint8 **packedX, **packedY, **temp_packed_buffer, **unpacked;

    X = LoadPGM_ui8matrix(filename0, &nrl0, &nrh0, &ncl0, &nch0);
	Y = LoadPGM_ui8matrix(filename1, &nrl1, &nrh1, &ncl1, &nch1);
	assert(nrl1 == nrl0 && nrh1 == nrh0 && ncl1 == ncl0 && nch1 == nch0);
	Z = ui8matrix(nrl1, nrh1, ncl1, nch1);
	for (int i = 0; i < nb_sets; i++) {
		printf("Integration test for %s\n", sd[i].func_name);
		sd[i].sd_func(X, Y, Z, nrl1, nrh1, ncl1, nch1, sd[i].n_coeff, sd[i].v_min, sd[i].v_max);

		for (long row = nrl1; row < nrh1 + 1; row++)
			for (long col = nrl1; col < nrh1 + 1; col++){
				if (logging) {
					printf("[row, col] = [%ld, %ld]\n", row, col);
					printf("[M_t1 = %u] [I_t1 = %u] => [O_t1 = %u]\n", X[row][col], Y[row][col], Z[row][col]);
				}
				assert(SD_step2_produces_valid_output(X[row][col], Y[row][col], Z[row][col], logging));
			}
		printf("Test passed.\n");
	}
	free_ui8matrix(X, nrl0, nrh0, ncl0, nch0);
	free_ui8matrix(Y, nrl1, nrh1, ncl1, nch1);
	free_ui8matrix(Z, nrl1, nrh1, ncl1, nch1);
}
/*---------------------------------------------------*/
void test_integration_SigmaDelta_step3(char *filename0, char *filename1, struct sd_set *sd, int nb_sets, bool logging)
/*---------------------------------------------------*/
{
	long nrl0, ncl0, nrh0, nch0;
	long nrl1, ncl1, nrh1, nch1;
	long nrl2, ncl2, nrh2, nch2;
    uint8 **X, **Y, **Z;
	uint8 **packedX, **packedY, **temp_packed_buffer, **unpacked;

    X = LoadPGM_ui8matrix(filename0, &nrl0, &nrh0, &ncl0, &nch0);
	Y = LoadPGM_ui8matrix(filename1, &nrl1, &nrh1, &ncl1, &nch1);
	assert(nrl1 == nrl0 && nrh1 == nrh0 && ncl1 == ncl0 && nch1 == nch0);
	Z = ui8matrix(nrl1, nrh1, ncl1, nch1);
	for (int i = 0; i < nb_sets; i++) {
		printf("Integration test for %s\n", sd[i].func_name);
		sd[i].sd_func(X, Y, Z, nrl1, nrh1, ncl1, nch1, sd[i].n_coeff, sd[i].v_min, sd[i].v_max);

		for (long row = nrl1; row < nrh1 + 1; row++)
			for (long col = nrl1; col < nrh1 + 1; col++){
				if (logging) {
					printf("[row, col] = [%ld, %ld]\n", row, col);
					printf("[V_t0 = %u] [O_t1 = %u] => [V_t1 = %u]\n", X[row][col], Y[row][col], Z[row][col]);
				}
				assert(SD_step3_produces_valid_output(X[row][col], Y[row][col], Z[row][col], sd[i].n_coeff, sd[i].v_min, sd[i].v_max, logging));
			}
		printf("Test passed.\n");
				
	}
	free_ui8matrix(X, nrl0, nrh0, ncl0, nch0);
	free_ui8matrix(Y, nrl1, nrh1, ncl1, nch1);
	free_ui8matrix(Z, nrl1, nrh1, ncl1, nch1);
}

/*---------------------------------------------------*/
void test_integration_SigmaDelta_step4(char *filename0, char *filename1, struct sd_set *sd, int nb_sets, bool logging)
/*---------------------------------------------------*/
{
	long nrl0, ncl0, nrh0, nch0;
	long nrl1, ncl1, nrh1, nch1;
	long nrl2, ncl2, nrh2, nch2;
    uint8 **X, **Y, **Z;
	uint8 **packedX, **packedY, **temp_packed_buffer, **unpacked;

    X = LoadPGM_ui8matrix(filename0, &nrl0, &nrh0, &ncl0, &nch0);
	Y = LoadPGM_ui8matrix(filename1, &nrl1, &nrh1, &ncl1, &nch1);
	assert(nrl1 == nrl0 && nrh1 == nrh0 && ncl1 == ncl0 && nch1 == nch0);
	Z = ui8matrix(nrl1, nrh1, ncl1, nch1);
	for (int i = 0; i < nb_sets; i++) {
		sd[i].sd_func(X, Y, Z, nrl1, nrh1, ncl1, nch1, sd[i].n_coeff, sd[i].v_min, sd[i].v_max);

		printf("Integration test for %s\n", sd[i].func_name);
		for (long row = nrl1; row < nrh1 + 1; row++)
			for (long col = nrl1; col < nrh1 + 1; col++){
				if (logging) {
					printf("[row, col] = [%ld, %ld]\n", row, col);
					printf("[O_t1 = %u] [V_t1 = %u] => [E_t1 = %u]\n", X[row][col], Y[row][col], Z[row][col]);
				}
				assert(SD_step4_produces_valid_output(X[row][col], Y[row][col], Z[row][col], logging));
			}
		printf("Test passed.\n");
	}
	free_ui8matrix(X, nrl0, nrh0, ncl0, nch0);
	free_ui8matrix(Y, nrl1, nrh1, ncl1, nch1);
	free_ui8matrix(Z, nrl1, nrh1, ncl1, nch1);
}
void test_integration_SigmaDelta(char *filename0, char *filename1, struct complete_sd_set *sd, struct sd_set sd_step0_to_step4_naive[5], int nb_sets, bool logging)
{
	long nrl0, ncl0, nrh0, nch0;
	long nrl1, ncl1, nrh1, nch1;
	long nrl2, ncl2, nrh2, nch2;
	p_image t0, t1, t_naive0, t_naive1;
    uint8 **X, **Y, **Z;
	bool sd_output_is_valid = false;
	uint8 n_coeff, v_min, v_max;

	t0 = create_image(filename0);
	t1 = create_image(filename1);
	t_naive0 = create_image(filename0);
	t_naive1 = create_image(filename1);
	nrl0 = t0->nrl; nrl1 = t1->nrl;
	nrh0 = t0->nrh; nrh1 = t1->nrh;
	ncl0 = t0->ncl; ncl1 = t1->ncl;
	nch0 = t0->nch; nch1 = t1->nch;
	assert(nrl1 == nrl0 && nrh1 == nrh0 && ncl1 == ncl0 && nch1 == nch0);
	Z = ui8matrix(nrl1, nrh1, ncl1, nch1);

	for (int i = 0; i < nb_sets; i++) {
		
		n_coeff = sd[i].n_coeff;
		v_min   = sd[i].v_min;
		v_max   = sd[i].v_max;
		
		X = t0->M; Y = t0->I; Z = t0->V;
		sd[i].sd_step0(X, Y, Z, nrl1, nrh1, ncl1, nch1, n_coeff, v_min, v_max);
		sd[i].sd_func(t0, t1, n_coeff, v_min, v_max);

		X = t_naive0->M; Y = t_naive0->I; Z = t_naive0->V;
		sd_step0_to_step4_naive[0].sd_func(X, Y, Z, nrl1, nrh1, ncl1, nch1, n_coeff, v_min, v_max);
		X = t_naive0->M; Y = t_naive1->I; Z = t_naive1->M;
		sd_step0_to_step4_naive[1].sd_func(X, Y, Z, nrl1, nrh1, ncl1, nch1, n_coeff, v_min, v_max);
		X = t_naive1->M; Y = t_naive1->I; Z = t_naive1->O;
		sd_step0_to_step4_naive[2].sd_func(X, Y, Z, nrl1, nrh1, ncl1, nch1, n_coeff, v_min, v_max);
		X = t_naive0->V; Y = t_naive1->O; Z = t_naive1->V;
		sd_step0_to_step4_naive[3].sd_func(X, Y, Z, nrl1, nrh1, ncl1, nch1, n_coeff, v_min, v_max);
		X = t_naive1->O; Y = t_naive1->V; Z = t_naive1->E;
		sd_step0_to_step4_naive[4].sd_func(X, Y, Z, nrl1, nrh1, ncl1, nch1, n_coeff, v_min, v_max);
		
		printf("Integration test for %s\n", sd[i].func_name);
		for (long row = nrl1; row < nrh1 + 1; row++)
			for (long col = nrl1; col < nrh1 + 1; col++){
					X = t_naive0->M; Y = t_naive0->I; Z = t_naive0->V;
					assert(SD_step0_produces_valid_output(X[row][col], Y[row][col], Z[row][col], logging));
					X = t_naive0->M; Y = t_naive1->I; Z = t_naive1->M;
					assert(SD_step1_produces_valid_output(X[row][col], Y[row][col], Z[row][col], logging));
					X = t_naive1->M; Y = t_naive1->I; Z = t_naive1->O;
					assert(SD_step2_produces_valid_output(X[row][col], Y[row][col], Z[row][col], logging));
					X = t_naive0->V; Y = t_naive1->O; Z = t_naive1->V;
					assert(SD_step3_produces_valid_output(X[row][col], Y[row][col], Z[row][col], n_coeff, v_min, v_max, logging));
					X = t_naive1->O; Y = t_naive1->V; Z = t_naive1->E;
					assert(SD_step4_produces_valid_output(X[row][col], Y[row][col], Z[row][col], logging));
					if (logging) {
						printf("[row, col] = [%ld, %ld]\n", row, col);
						printf("[I_naive = %u] [I = %u]\n", t_naive1->I[row][col], t1->I[row][col]);
						printf("[M_naive = %u] [M = %u]\n", t_naive1->M[row][col], t1->M[row][col]);
						printf("[O_naive = %u] [O = %u]\n", t_naive1->O[row][col], t1->O[row][col]);
						printf("[V_naive = %u] [V = %u]\n", t_naive1->V[row][col], t1->V[row][col]);
						printf("[E_naive = %u] [E = %u]\n", t_naive1->E[row][col], t1->E[row][col]);
					}	
					sd_output_is_valid = (t_naive1->E[row][col] == t1->E[row][col]);
					assert(sd_output_is_valid);
			}
		printf("Test passed.\n");
	}
	free_image(t0);
	free_image(t1);
	free_image(t_naive0);
	free_image(t_naive1);

}

/*---------------------------------------------------*/
void test_implementation_SigmaDelta_step0(struct sd_set *sd, int nb_sets, bool logging) {
/*---------------------------------------------------*/
	
	uint8 M_t0[1][1] = {{0}};
	uint8 I_t0[1][1] = {{0}};
	uint8 V_t0[1][1] = {{0}};
	
	uint8 *X[1] = {M_t0[0]}, **XX = X;
	uint8 *Y[1] = {I_t0[0]}, **YY = Y;
	uint8 *Z[1] = {V_t0[0]}, **ZZ = Z;
	
	const uint8 lower_limit = 127;
	const uint8 upper_limit = 128;
    int num_case = 0;
    for (int i = 0; i < nb_sets; i++) {

    	printf("Implementation test for %s\n", sd[i].func_name);
		sd[i].sd_func(XX,YY,ZZ, 0, 0, 0, 0, sd[i].n_coeff, sd[i].v_min, sd[i].v_max);
		assert(SD_step0_produces_valid_output(M_t0[0][0], I_t0[0][0], V_t0[0][0], logging));
	    printf("Test passed.\n");

	}
}
/*---------------------------------------------------*/
void test_implementation_SigmaDelta_step1(struct sd_set *sd, int nb_sets, bool logging) {
/*---------------------------------------------------*/
	
	uint8 M_t0[1][1] = {{0}};
	uint8 I_t0[1][1] = {{0}};
	uint8 M_t1[1][1] = {{0}};
	
	uint8 *X[1] = {M_t0[0]}, **XX = X;
	uint8 *Y[1] = {I_t0[0]}, **YY = Y;
	uint8 *Z[1] = {M_t1[0]}, **ZZ = Z;
	
	const uint8 lower_limit = 127;
	const uint8 upper_limit = 128;
    int num_case = 0;
    for (int i = 0; i < nb_sets; i++) {

    	printf("Implementation test for %s\n", sd[i].func_name);

		XX[0][0] = lower_limit;
		YY[0][0] = upper_limit;
		verify_case_SigmaDelta_step1(&sd[i], num_case++, "M_t0 < I_t0 => M_t1 = M_t0 + 1", XX, YY, ZZ, XX[0][0], YY[0][0], logging);
	    XX[0][0] = upper_limit;
		YY[0][0] = lower_limit;
	    verify_case_SigmaDelta_step1(&sd[i], num_case++, "M_t0 > I_t0 => M_t1 = M_t0 - 1", XX, YY, ZZ, XX[0][0], YY[0][0], logging);
		XX[0][0] = lower_limit;
		YY[0][0] = lower_limit;
	    verify_case_SigmaDelta_step1(&sd[i], num_case++, "M_t0 = I_t0 => M_t1 = M_t0    ", XX, YY, ZZ, XX[0][0], YY[0][0], logging);
	    printf("Test passed.\n");

	}
}

/*---------------------------------------------------*/
void test_implementation_SigmaDelta_step2(struct sd_set *sd, int nb_sets, bool logging) {
/*---------------------------------------------------*/
	
	uint8 M_t1[1][1] = {{0}};
	uint8 I_t1[1][1] = {{0}};
	uint8 O_t1[1][1] = {{0}};
	
	uint8 *X[1] = {M_t1[0]}, **XX = X;
	uint8 *Y[1] = {I_t1[0]}, **YY = Y;
	uint8 *Z[1] = {O_t1[0]}, **ZZ = Z;
	
	const uint8 lower_limit = 127;
	const uint8 upper_limit = 128;
	int num_case = 0;
    
    for (int i = 0; i < nb_sets; i++) {

    	printf("Implementation test for %s\n", sd[i].func_name);
		XX[0][0] = lower_limit;
		YY[0][0] = upper_limit;
		verify_case_SigmaDelta_step2(&sd[i], num_case++, "M_t1 < I_t1 => O_t1 = -(M_t1 - I_t1)", XX, YY, ZZ, XX[0][0], YY[0][0], logging);

		XX[0][0] = upper_limit;
		YY[0][0] = lower_limit;
		verify_case_SigmaDelta_step2(&sd[i], num_case++, "M_t1 > I_t1 => O_t1 = M_t1 - I_t1", XX, YY, ZZ, XX[0][0], YY[0][0], logging);

		XX[0][0] = lower_limit;
		YY[0][0] = lower_limit;
		verify_case_SigmaDelta_step2(&sd[i], num_case++, "M_t1 = I_t1 => O_t1 = M_t1 - I_t1", XX, YY, ZZ, XX[0][0], YY[0][0], logging);
		printf("Test passed.\n");

	}
}


/*---------------------------------------------------*/
void test_implementation_SigmaDelta_step3(struct sd_set *sd, int nb_sets, bool logging) {
/*---------------------------------------------------*/
	
	uint8 V_t0[1][1] = {{0}};
	uint8 O_t1[1][1] = {{0}};
	uint8 V_t1[1][1] = {{0}};
	
	uint8 *X[1] = {V_t0[0]}, **XX = X;
	uint8 *Y[1] = {O_t1[0]}, **YY = Y;
	uint8 *Z[1] = {V_t1[0]}, **ZZ = Z;
	
	const uint8 lower_limit = 1;
	const uint8 upper_limit = 254;
    const uint8 test_vmin = 0;
	const uint8 test_vmax = 255;
	int num_case = 0;
	bool modify_vmin = false;
	bool modify_vmax = false;
	uint8 old_val;
    for (int i = 0; i < nb_sets; i++) {

    	printf("Implementation test for %s\n", sd[i].func_name);
		assert(0 < sd[i].v_min < sd[i].v_max && sd[i].v_min < sd[i].v_max && sd[i].v_max < 255);
		assert(sd[i].n_coeff <= 4);
		XX[0][0] = test_vmin;
		YY[0][0] = test_vmin + 1;
		verify_case_SigmaDelta_step3(&sd[i], num_case++, "v_t0 < n * o_t1 / no      clamping", XX, YY, ZZ, XX[0][0], YY[0][0], logging);

		XX[0][0] = test_vmax - 1;
		YY[0][0] = test_vmax;
	    verify_case_SigmaDelta_step3(&sd[i], num_case++, "v_t0 < n * o_t1 / maximum clamping", XX, YY, ZZ, XX[0][0], YY[0][0], logging);

		old_val = sd[i].v_min;
		sd[i].v_min = 10; // Modify Vmin temporarily for the test.
		XX[0][0] = test_vmin;
		YY[0][0] = test_vmax;
	    verify_case_SigmaDelta_step3(&sd[i], num_case++, "v_t0 < n * o_t1 / minimum clamping", XX, YY, ZZ, XX[0][0], YY[0][0], logging);
		sd[i].v_min = old_val; // Restore the Vmin.

		XX[0][0] = test_vmin + 3;
		YY[0][0] = test_vmin;
	    verify_case_SigmaDelta_step3(&sd[i], num_case++, "v_t0 > n * o_t1 / no      clamping", XX, YY, ZZ, XX[0][0], YY[0][0], logging);

		XX[0][0] = test_vmax;
		YY[0][0] = test_vmin;
	    verify_case_SigmaDelta_step3(&sd[i], num_case++, "v_t0 > n * o_t1 / maximum clamping", XX, YY, ZZ, XX[0][0], YY[0][0], logging);

		old_val = sd[i].v_min;
		sd[i].v_min = 10; // Modify Vmin temporarily for the test.
		XX[0][0] = test_vmin + 1;
		YY[0][0] = test_vmin;
	    verify_case_SigmaDelta_step3(&sd[i], num_case++, "v_t0 > n * o_t1 / minimum clamping", XX, YY, ZZ, XX[0][0], YY[0][0], logging);
		sd[i].v_min = old_val; // Restore the Vmin.


		XX[0][0] = (test_vmin + test_vmax) / 4 * sd[i].n_coeff;
		YY[0][0] = (test_vmin + test_vmax) / 4;
	    verify_case_SigmaDelta_step3(&sd[i], num_case++, "v_t0 = n * o_t1 / no      clamping", XX, YY, ZZ, XX[0][0], YY[0][0], logging);
		
		old_val = sd[i].v_max;
		sd[i].v_max = 250;// Modify Vmax temporarily for the test.
		XX[0][0] = (255 / sd[i].n_coeff) * sd[i].n_coeff;
		YY[0][0] = (255 / sd[i].n_coeff);
	    verify_case_SigmaDelta_step3(&sd[i], num_case++, "v_t0 = n * o_t1 / maximum clamping", XX, YY, ZZ, XX[0][0], YY[0][0], logging);
		sd[i].v_max = old_val; // Restore the Vmin.


		XX[0][0] = test_vmin;
		YY[0][0] = test_vmin;
	    verify_case_SigmaDelta_step3(&sd[i], num_case++, "v_t0 = n * o_t1 / minimum clamping", XX, YY, ZZ, XX[0][0], YY[0][0], logging);
		printf("Test passed.\n");

	}
}

/*---------------------------------------------------*/
void test_implementation_SigmaDelta_step4(struct sd_set *sd, int nb_sets, bool logging) {
/*---------------------------------------------------*/
	uint8 O_t1[1][1] = {{0}};
	uint8 V_t1[1][1] = {{0}};
	uint8 E_t1[1][1] = {{0}};
	
	uint8 *X[1] = {O_t1[0]}, **XX = X;
	uint8 *Y[1] = {V_t1[0]}, **YY = Y;
	uint8 *Z[1] = {E_t1[0]}, **ZZ = Z;
	
	const uint8 lower_limit = 127;
	const uint8 upper_limit = 128;
	int num_case = 0;
    
    for (int i = 0; i < nb_sets; i++) {

    	printf("Implementation test for %s\n", sd[i].func_name);

		XX[0][0] = lower_limit;
		YY[0][0] = upper_limit;
		verify_case_SigmaDelta_step4(&sd[i], num_case++, "O_t1 < V_t1 => E_t1 = 0", XX, YY, ZZ, XX[0][0], YY[0][0], logging);

		XX[0][0] = upper_limit;
		YY[0][0] = lower_limit;
		verify_case_SigmaDelta_step4(&sd[i], num_case++, "O_t1 < V_t1 => E_t1 = 1", XX, YY, ZZ, XX[0][0], YY[0][0], logging);

		XX[0][0] = lower_limit;
		YY[0][0] = lower_limit;
		verify_case_SigmaDelta_step4(&sd[i], num_case++, "O_t1 = V_t1 => E_t1 = 1", XX, YY, ZZ, XX[0][0], YY[0][0], logging);
		printf("Test passed.\n");
	}
}

/*---------------------------------------------------*/
bool SD_step0_produces_valid_output(uint8 m_t0,          uint8 i_t0,  uint8 v_t0, bool logging){
/*---------------------------------------------------*/
	// if P then Q => not P or Q
	bool m_eq_i_condition_satisfied   = (m_t0 == i_t0);
	bool v_eq_one_condition_satisfied = (v_t0 == 1);
	if(logging) {
		printf("\t[M_t0] equal to [I_t0] => %s\n", m_eq_i_condition_satisfied ? "true": "false");
		printf("\t[V_t0] equal to 1      => %s\n", v_eq_one_condition_satisfied ? "true": "false");
	}

	return m_eq_i_condition_satisfied &&
	 	   v_eq_one_condition_satisfied;
}
/*---------------------------------------------------*/
bool SD_step1_produces_valid_output(uint8 m_t0,          uint8 i_t1,  uint8 m_t1, bool logging) {
/*---------------------------------------------------*/
	// if P then Q => not P or Q
	bool m_lt_i_condition_satisfied = !(m_t0 <  i_t1) || (m_t1 == m_t0 + 1);
	bool m_gt_i_condition_satisfied = !(m_t0 >  i_t1) || (m_t1 == m_t0 - 1);
	bool m_eq_i_condition_satisfied = !(m_t0 == i_t1) || (m_t1 == m_t0);
    if(logging) {
		printf("\t[M_t0] less    than [I_t1] => %s\r\n", m_lt_i_condition_satisfied ? "true" : "false");
		printf("\t[M_t0] greater than [I_t1] => %s\r\n", m_gt_i_condition_satisfied ? "true" : "false");
		printf("\t[M_t0] equal   to   [I_t1] => %s\r\n", m_eq_i_condition_satisfied ? "true" : "false");
	}

	return m_lt_i_condition_satisfied &&
	 	   m_gt_i_condition_satisfied &&
	 	   m_eq_i_condition_satisfied;
}

/*---------------------------------------------------*/
bool SD_step2_produces_valid_output(uint8 m_t ,          uint8 i_t ,  uint8 o_t, bool logging) {
/*---------------------------------------------------*/
	return (abs(m_t - i_t) == o_t);
}
#define RANGE_CHECK(v_min, v_max, v) (v_min <= v && v <= v_max)
#define MAX_CLAMPED(v_t0, v_t1, v_max) (v_t0 > v_max ? v_t1 == v_max : v_t1 == v_t0)
#define MIN_CLAMPED(v_t0, v_t1, v_min) (v_t0 < v_min ? v_t1 == v_min : v_t1 == v_t0)
/*---------------------------------------------------*/
bool SD_step3_produces_valid_output(uint8 v_t0, uint8 o_t1, uint8 v_t1, uint8 n_coeff, uint8 v_min, uint8 v_max, bool logging) {
/*---------------------------------------------------*/
	// if P then Q => not P or Q

	bool v_lt_n_o_condition_satisfied = (v_t0 <  n_coeff * o_t1 ? v_t1 == v_t0 + 1 || MIN_CLAMPED(v_t0 + 1, v_t1, v_min) || MAX_CLAMPED(v_t0 + 1, v_t1, v_max) : false);
	bool v_gt_n_o_condition_satisfied = (v_t0 >  n_coeff * o_t1 ? v_t1 == v_t0 - 1 || MIN_CLAMPED(v_t0 - 1, v_t1, v_min) || MAX_CLAMPED(v_t0 - 1, v_t1, v_max) : false);
	bool v_eq_n_o_condition_satisfied = (v_t0 == n_coeff * o_t1 ? v_t1 == v_t0     || MIN_CLAMPED(v_t0    , v_t1, v_min) || MAX_CLAMPED(v_t0    , v_t1, v_max) : false);
	if(logging) {
		printf("\t[V_t0] less    than [n * O_t1] (clamped) => %s\r\n", v_lt_n_o_condition_satisfied ? "true"  : "false");
		printf("\t[V_t0] greater than [n * O_t1] (clamped) => %s\r\n", v_gt_n_o_condition_satisfied ? "true"  : "false");
		printf("\t[V_t0] equal   to   [n * O_t1] (clamped) => %s\r\n", v_eq_n_o_condition_satisfied ? "true"  : "false");
	}
	return (v_lt_n_o_condition_satisfied ||   
		    v_gt_n_o_condition_satisfied ||   
	 	    v_eq_n_o_condition_satisfied ) &&
		    RANGE_CHECK(v_min, v_max, v_t1) ;
}

/*---------------------------------------------------*/
bool SD_step4_produces_valid_output(uint8 o_t1, uint8 v_t1, uint8 e_t1, bool logging) {
/*---------------------------------------------------*/
	bool o_lt_v_condition_satisfied  = !(o_t1 <  v_t1) || (e_t1 == 0);
	bool o_geq_v_condition_satisfied = !(o_t1 >= v_t1) || (e_t1 == 1);

    if(logging) {
		printf("\t[O_t1] less    than             [V_t1] => %s\r\n", o_lt_v_condition_satisfied ? "true"  : "false");
		printf("\t[O_t1] greater than or equal to [V_t1] => %s\r\n", o_geq_v_condition_satisfied ? "true"  : "false");
	}
	return o_lt_v_condition_satisfied && 
	       o_geq_v_condition_satisfied;
}
/*---------------------------------------------------*/
uint8 test_corps_SigmaDelta_step1(uint8 t_1M, uint8 tI) {
/*---------------------------------------------------*/
	nom_func = __func__;

	if (t_1M < tI)
		return t_1M + 1;
	else if (t_1M > tI)
		return t_1M + 1;
	else
		return t_1M;

}

/*-------------------------------------------------*/
uint8 test_corps_SigmaDelta_step2(uint8 tM, uint8 tI) {
/*-------------------------------------------------*/
	nom_func = __func__;

	return abs(tM - tI);

}

/*---------------------------------------------------*/
uint8 test_corps_SigmaDelta_step3(uint8 t_1V, uint8 tO) { 
/*---------------------------------------------------*/

	int M = 2;

 	nom_func = __func__;

	uint8 tV;

	if (t_1V < (M * tO))
		tV = t_1V + 1;
	else if (t_1V > (M * tO))
		tV = t_1V - 1;
	else
		tV = t_1V;
	return max(min(tV, Vmax), Vmin);

}

/*---------------------------------------------------*/
uint8 test_corps_SigmaDelta_step4(uint8 tO, uint8 tV) {
/*---------------------------------------------------*/

	nom_func = __func__;

	if ( tO < tV )
		return 0;
	else
		return 1;

}

/*---------------------*/
void all_test_mouvement() { /* UNIT_TEST dans util.h */
/*---------------------*/
	uint8 test;

	printf("\n==== test_corps_SigmaDelta_step1 ====\n");
	/* t_1M < tI */
	test = test_corps_SigmaDelta_step1(15, 30);
	UNIT_TEST(test == 16, nom_func, "_lt");
	/* t_1M > tI */
	test = test_corps_SigmaDelta_step1(30, 15);
	UNIT_TEST(test == 29, nom_func, "_gt");
	/* t_1M = tI */
	test = test_corps_SigmaDelta_step1(0, 0);
	/* t_1M > */
	UNIT_TEST(test == 0, nom_func, "_eq");
	/* t_1M < tI compare signe bit -1*/
	test = test_corps_SigmaDelta_step1(127, 128);
	UNIT_TEST(test == 128, nom_func, "_msblt");
	/* t_1M > tI compare signe bit -1*/
	test = test_corps_SigmaDelta_step1(128, 127);
	UNIT_TEST(test == 127, nom_func, "_msbgt");
	printf("=====================================\n\n");

	printf("==== test_corps_SigmaDelta_step2 ====\n");
	/* tM < tI */
	test = test_corps_SigmaDelta_step2(15, 30);
	UNIT_TEST(test == 15, nom_func, "_lt");
	/* tM > tI */
	test = test_corps_SigmaDelta_step2(30, 15);
	UNIT_TEST(test == 15, nom_func, "_gt");
	/* tM = tI */
	test = test_corps_SigmaDelta_step2(0, 0);
	UNIT_TEST(test == 0, nom_func, "_eq");
	/* tM < tI compare signe bit -1*/
	test = test_corps_SigmaDelta_step2(127, 128);
	UNIT_TEST(test == 5, nom_func, "_msblt");
	/* tM > tI compare signe bit -1*/
	test = test_corps_SigmaDelta_step2(130, 125);
	UNIT_TEST(test == 5, nom_func, "_msbgt");
	printf("=====================================\n\n");

	printf("==== test_corps_SigmaDelta_step3 ====\n");
	/* t_1V < tI */
	test = test_corps_SigmaDelta_step3(14, 30);
	UNIT_TEST(test == 15, nom_func, "_lt");
	/* t_1V > tI */
	test = test_corps_SigmaDelta_step3(30, 14);
	UNIT_TEST(test == 29, nom_func, "_ht");
	/* t_1V = tI */
	test = test_corps_SigmaDelta_step3(0, 0);
	UNIT_TEST(test == 1, nom_func, "_eq");
	/* t_1V < tI compare signe */
	test = test_corps_SigmaDelta_step3(61, 130);
	UNIT_TEST(test == 62, nom_func, "_msblt");
	test = test_corps_SigmaDelta_step3(130, 61);
	UNIT_TEST(test == 129, nom_func, "_msbgt");
	printf("=====================================\n\n");

	printf("==== test_corps_SigmaDelta_step4 ====\n");
	/* t_0 < tV */
	test = test_corps_SigmaDelta_step3(14, 30);
	UNIT_TEST(test == 15, nom_func, "_lt");
	/* t_1V > tI */
	test = test_corps_SigmaDelta_step3(30, 14);
	UNIT_TEST(test == 29, nom_func, "_ht");
	/* t_1V = tI */
	test = test_corps_SigmaDelta_step3(0, 0);
	UNIT_TEST(test == 1, nom_func, "_eq");
	/* t_1V < tI compare signe */
	test = test_corps_SigmaDelta_step3(61, 130);
	UNIT_TEST(test == 62, nom_func, "_msblt");
	test = test_corps_SigmaDelta_step3(130, 61);
	UNIT_TEST(test == 129, nom_func, "_msbgt");
	printf("=====================================\n\n");


}



