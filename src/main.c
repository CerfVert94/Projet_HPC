/* -------------- */
/* --- main.c --- */
/* -------------- */

#include <stdio.h>
#include <stdlib.h>
#include <emmintrin.h>
#include <tmmintrin.h>

#include "nrdef.h"
#include "nrutil.h"

#include "vnrdef.h"
#include "vnrutil.h"

#include "util.h"

#include "mymacro.h"
#include "simd_macro.h"

#include "img.h"
#include "img_SIMD.h"
#include "mouvement.h"
#include "mouvement_SIMD.h"
#include "morpho.h"
#include "morpho_SIMD.h"
#include "test_mouvement.h"

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

void test_val() {
    int card = card_vuint8 ();

    int v1, v2;
    int m1, m2;

    long rl, rh, cl, ch;

    uint8** tmp;
    tmp = LoadPGM_ui8matrix("../car3/car_3000.pgm", &rl, &rh, &cl, &ch);

    s2v1D(rl, rh, card, &v1, &v2);
    v2m1D(v1, v2, card, &m1, &m2);
    printf("%ld %ld \n%d %d\n%d %d\n", rl, rh, v1, v2, m1, m2);
    vuint8 x;

    vuint8** vtmp;
    vtmp = vui8matrix(rl, rh, v1, v2);

    int i, j, k, l;

    for (i = rl; i < rh; i++) {
        for (j = v1; j <= v2; j++) {
            //x = init_vuint8(11);
            l = j*card;
            x = init_vuint8_all(tmp[i][l], tmp[i][l+1], tmp[i][l+2], tmp[i][l+3], tmp[i][l+4], tmp[i][l+5], tmp[i][l+6], tmp[i][l+7], tmp[i][l+8], tmp[i][l+9], tmp[i][l+10], tmp[i][l+11], tmp[i][l+12], tmp[i][l+13], tmp[i][l+14], tmp[i][l+15]);
            _mm_store_si128(&vtmp[i][j], x);
        }
    }

    display_ui8vector((uint8*) tmp[0], rl, rh, "%4d", "TMP"); puts("");
    display_vui8vector(vtmp[0], v1, v2, "%4d", "vTMP"); puts("");
    
    free_ui8matrix(tmp, rl, rh, cl, ch);
    free_vui8matrix(vtmp, rl, rh, v1, v2);
}

void test_abs() {

    vuint8 A, B, C, D, CMP;

    A = init_vuint8(-10);
    B = init_vuint8(-50);
    //C = init_vuint8_all(4, 1, 4, 1, 4, 1, 4, 1, 4, 1, 4, 1, 4, 1, 4, 1);
    CMP = init_vuint8(127);

    D = vec_subabs(B,A);

    display_vuint8(D, "%4d", NULL);

}

void test_store_si128() {

    vuint8 A;

    vuint8 B = init_vuint8(5);
    _mm_store_si128((vuint8*) &A, B);

    display_vuint8(A, "%4d", NULL);
}

void test_gt_eq() {

    vuint8 A, B, C1, C2, D, E, F, G, CMP;

    vuint8 AB;

    A = init_vuint8(254);
    B = init_vuint8(255);

    /*E = init_vuint8(6);
    F = init_vuint8(3);
    G = init_vuint8(19);

    AB = F;

    for (int i = 1; i < N; i++) {
        AB = _mm_add_epi8(AB, F);
    }
    //AB = _mm_mulhi_epi16(AB, vN);
    display_vuint8(AB, "%4d", NULL);
    */

    // if (A < B) D = 6; else d = 3;

    CMP = init_vuint8(128);

    //vec_cmplt(A, B, C1, CMP);

    A = _mm_sub_epi8(A, CMP);
    B = _mm_sub_epi8(B, CMP);
    C1 = _mm_cmplt_epi8(A, B);
    A = _mm_add_epi8(A, CMP);
    B = _mm_add_epi8(B, CMP);
    //C2  = _mm_cmpeq_epi8(A, B);

    display_vuint8(C1, "%4d", NULL);
    //display_vuint8(C2, "%4d", NULL);
    display_vuint8(A, "%4d", NULL);
    display_vuint8(B, "%4d", NULL);

    /*D = _mm_or_si128(_mm_and_si128(C1, E), _mm_andnot_si128(C1, F));
    display_vuint8(D, "%4d", NULL);
    D = _mm_or_si128(_mm_and_si128(C2, G), _mm_andnot_si128(C2, D));
    display_vuint8(D, "%4d", NULL);*/

}

void test_if_then_else() {

    vuint8 A, B, C, D, E, CMP;

    vuint8 GT, LT;

    vsint8 TEST = init_vsint8(-1);
    //display_vsint8(TEST, "%4d", NULL);


    GT = init_vuint8(1);
    LT = init_vuint8(0);
    A = init_vuint8(255);
    B = init_vuint8(127);

    CMP = init_vuint8(127);

    A = _mm_sub_epi8(A, CMP);
    B = _mm_sub_epi8(B, CMP);
    C = _mm_cmpeq_epi8(A, B);
    A = _mm_add_epi8(A, CMP);
    B = _mm_add_epi8(B, CMP);

    display_vuint8(C, "%4d", NULL);
    display_vuint8(A, "%4d", NULL);
    display_vuint8(B, "%4d", NULL);

    D = _mm_or_si128(_mm_and_si128(C, LT), _mm_andnot_si128(C, GT));
    display_vuint8(D, "%4d", NULL);

}

void test_signed() {

    vuint8 A, C, D, E, CMP;

    vsint8 B;

    A = init_vuint8(127);

    B = init_vsint8(-1);

    C = init_vuint8(1);

    D = _mm_add_epi8(A, B);

    display_vuint8(D, "%4d", NULL);

    E = _mm_add_epi8(A, C);

    display_vuint8(E, "%4d", NULL);


}

void test_macro_shift() {

    vuint8 A = init_vuint8(4);
    vuint8 B = init_vuint8(2);
    vuint8 C = init_vuint8(3);

    vuint8 D = init_vuint8(5);
    vuint8 E = init_vuint8(6);
    vuint8 F = init_vuint8(5);

    vuint8 G = init_vuint8(10);
    vuint8 H = init_vuint8(56);
    vuint8 I = init_vuint8(100);

    vuint8 X, Y, Z;
    vuint8 OM;
    X = vector_and3(A,B,C);
    display_vuint8(X, "%4d", NULL);
    printf("\n");
    Y = vector_and3(D,E,F);
    display_vuint8(Y, "%4d", NULL);
    printf("\n");
    Z = vector_and3(G,H,I);
    display_vuint8(Z, "%4d", NULL);
    printf("\n");
    printf("\n\n");

    OM = vec_right1(X, Y);
    display_vuint8(OM, "%4d", NULL);
    printf("\n");
    OM = vec_left1(Y, Z);
    display_vuint8(OM, "%4d", NULL);
    printf("\n");
    OM = vector_and3(vec_right1(X, Y),Y,vec_left1(Y, Z));
    display_vuint8(OM, "%4d", NULL);
    printf("\n");
}

void test_cast() {

    uint8* X = (uint8*)malloc(32*sizeof(uint8));
    vuint8 vX;
    vX = init_vuint8_all(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16);

    X = (uint8*)&vX;

    for(int i = 0; i < 32; i++){
        printf("%d ", X[i]);
    }
    puts("");

}

// -----------
int main(void)
// -----------
{
    //info();
    //all_test_mouvement();
    //test_mouvement();

    //test_step0_SIMD();
    test2_SIMD_img();
    //test_abs();
    //test_val();
    //test_store_si128();
    //test_gt_eq();
    //test_if_then_else();
    //test_signed();
    //test_macro_shift();
    //test_functions_morpho_SIMD();
    //test_cast();
    return 0;    
}