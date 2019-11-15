/* -------------- */
/* --- main.c --- */
/* -------------- */

#include <stdio.h>
#include <stdlib.h>


#include "nrdef.h"
#include "nrutil.h"

#include "vnrdef.h"
#include "vnrutil.h"


#include "img.h"
#include "mouvement.h"
#include "mymacro.h"
//#include "test_simd1.h"
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

// -----------
int main(void)
// -----------
{
    //info();
    //all_test_mouvement();

    int card = card_vuint8 ();

    int v1, v2;
    int m1, m2;

    long rl, rh, cl, ch;

    uint8** tmp;
    tmp = LoadPGM_ui8matrix("../car3/car_3000.pgm", &rl, &rh, &cl, &ch);

    s2v1D(rl, rh, card, &v1, &v2);
    v2m1D(v1, v2, card, &m1, &m2);
    //printf("%ld %ld \n%d %d\n%d %d\n", rl, rh, v1, v2, m1, m2);

    vuint8 x;

    vuint8** vtmp;
    vtmp = vui8matrix(rl, rh, v1, v2);

    int i, j, k, l;

    for (i = rl; i < rh; i++) {
        for (j = v1; j <= v2; j++) {
            //x = init_vuint8(11);
            l = j*16;
            x = init_vuint8_all(tmp[i][j], tmp[i][l+1], tmp[i][l+2], tmp[i][l+3], tmp[i][l+4], tmp[i][l+5], tmp[i][l+6], tmp[i][l+7], tmp[i][l+8], tmp[i][l+9], tmp[i][l+10], tmp[i][l+11], tmp[i][l+12], tmp[i][l+13], tmp[i][l+14], tmp[i][l+15]);
            _mm_store_si128(&vtmp[i][j], x);
        }
    }

    display_ui8vector((uint8*) vtmp[65], rl, rh, "%4d", "TMP"); puts("");
    display_vui8vector(vtmp[65], v1, v2, "%4d", "vTMP"); puts("");
    
    free_ui8matrix(tmp, rl, rh, cl, ch);
    free_vui8matrix(vtmp, rl, rh, v1, v2);


    return 0;    
}