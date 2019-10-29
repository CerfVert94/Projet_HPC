 /* --------------- */
/* --- simd1.c --- */
/* --------------- */

#include <stdio.h>
#include <stdlib.h>

#include "nrdef.h"
#include "nrutil.h"

#include "vnrdef.h"
#include "vnrutil.h"

#include "mutil.h"

#include "mymacro.h"
#include "simd_macro.h"
#include "simd1.h"

// ------------------
void test_macro(void)
// ------------------
{
    vfloat32 a, b, c, d, e;
    vfloat32 au, cu; // unaligned vector
    vfloat32 a3, a5; // add3 add5
    puts("------------------");
    puts("--- test_macro ---");
    puts("------------------"); puts("");
    
    puts("vec_left et vec_right");
    
    // setr simule le fonctionnement d'un load avec permutation des blocs grace au "r" = reverse
    a = _mm_setr_ps(1,   2,  3,  4);
    b = _mm_setr_ps(5,   6,  7,  8);
    c = _mm_setr_ps(9,  10, 11, 12);
    d = _mm_setr_ps(13, 14, 15, 16);
    e = _mm_setr_ps(17, 18, 19, 20);
    
    display_vfloat32(a, "%4.0f", "a "); puts("");
    display_vfloat32(b, "%4.0f", "b "); puts("");
    display_vfloat32(c, "%4.0f", "c "); puts("\n");
    
    au = vec_left1(a, b);
    cu = vec_right1(b, c);
    display_vfloat32(au, "%4.0f", "au"); puts("");
    display_vfloat32(cu, "%4.0f", "cu"); puts("\n");
    
    puts("vec_add3 et vec_add5");
    a3 = vec_add3(a, b, c);
    a5 = vec_add5(a, b, c, d, e);
    
    display_vfloat32(a3, "%4.0f", "a3"); puts("");
    display_vfloat32(a5, "%4.0f", "a5"); puts("\n");
}
// -------------------------------------------------------------------
void add_vf32vector(vfloat32 *vX1, vfloat32 *vX2, int n, vfloat32 *vY)
// -------------------------------------------------------------------
{
    vfloat32 x1, x2, y;

    for(int i=0; i<n; i++) {
        
        x1 = _mm_load_ps((float32*) &vX1[i]);
        x2 = _mm_load_ps((float32*) &vX2[i]);
        
        y = _mm_add_ps(x1, x2);
        
        DEBUG(display_vfloat32(x1, "%4.0f", "x1 =")); DEBUG(puts(""));
        DEBUG(display_vfloat32(x2, "%4.0f", "x2 =")); DEBUG(puts(""));
        DEBUG(display_vfloat32(y,  "%4.0f", "y  =")); DEBUG(puts(""));
        
        _mm_store_ps((float*) &vY[i], y);
        
        DEBUG(puts("-------------------"));
    }
}
// ---------------------------------------------------------
vfloat32 dot_vf32vector(vfloat32 *vX1, vfloat32 *vX2, int n)
// ---------------------------------------------------------
{
    vfloat32 z;
    vfloat32 x1, x2, p, s;

    s = _mm_set_ps(0,0,0,0);

    
    
    return s; // attention il faut retourner un registre SIMD et non un scalaire

}
// ----------------------------------------------------
void sum3_vf32vector(vfloat32 *vX, int n, vfloat32 *vY)
// ----------------------------------------------------
{
}
// ----------------------------------------------------
void sum5_vf32vector(vfloat32 *vX, int n, vfloat32 *vY)
// ----------------------------------------------------
{
}
// ----------------------------------------------------
void min3_vf32vector(vfloat32 *vX, int n, vfloat32 *vY)
// ----------------------------------------------------
{
}
// ----------------------------------------------------
void min5_vf32vector(vfloat32 *vX, int n, vfloat32 *vY)
// ----------------------------------------------------
{
}
// -------------------------------------------------------------
void positive_add3_vf32vector(vfloat32 *vX, int n, vfloat32 *vY)
// -------------------------------------------------------------
{
}
// -------------------------------------------------------------
void positive_avg3_vf32vector(vfloat32 *vX, int n, vfloat32 *vY)
// -------------------------------------------------------------
{
}
// --------------------------------------------------
vfloat32 positive_avg_vf32vector(vfloat32 *vX, int n)
// --------------------------------------------------
{
    vfloat32 avg3 = _mm_set1_ps(0.0f);
    
    return avg3;

}