/* -------------------- */
/* --- simd_macro.h --- */
/* -------------------- */


#ifndef __SIMD_MACRO_H__
#define __SIMD_MACRO_H__

#pragma message("  include  simd_macro.h")

// -----------------------
// melanges / permutations
// -----------------------

#define vec_border(a,b) b
#define vec_middle(a,b) a

#define vec_left1(v0, v1) _mm_shuffle_ps(v0, _mm_shuffle_ps(v0, v1, _MM_SHUFFLE(1,0,3,2)), _MM_SHUFFLE(2,1,2,1))
#define vec_left2(v0, v1) _mm_shuffle_ps(v0, v1, _MM_SHUFFLE(1,0,3,2))
#define vec_left3(v0, v1) _mm_shuffle_ps(_mm_shuffle_ps(v0, v1, _MM_SHUFFLE(2,0,3,2)), v1, _MM_SHUFFLE(2,1,2,1))
#define vec_left4(v0, v1) v1

#define vec_right1(v1, v2) _mm_shuffle_ps(_mm_shuffle_ps(v2, v1, _MM_SHUFFLE(3,2,1,0)), v2, _MM_SHUFFLE(2,1,0,3))
#define vec_right2(v1, v2) _mm_shuffle_ps(v1, v2, _MM_SHUFFLE(1,0,3,2))
#define vec_right3(v1, v2) _mm_shuffle_ps(v1, _mm_shuffle_ps(v2, v1, _MM_SHUFFLE(3,2,1,0)), _MM_SHUFFLE(0,3,2,1))
#define vec_right4(v1, v2) v1 

// -------
//   MIN
// -------
#define vMIN2(x0,x1) x0  
#define vMIN3(x0,x1,x2) x0       
#define vMIN4(x0,x1,x2,x3) x0  
#define vMIN5(x0,x1,x2,x3,x4) x0

// -------
// calculs
// -------

#define vec_add3(x0, x1, x2) x0
#define vec_add5(x0, x1, x2, x3, x4) x0

#define vec_min3(x0, x1, x2) x0
#define vec_min5(x0, x1, x2, x3, x4) x0

// division neutralisee pour verifier la somme
#define vec_div3(x) x
#define vec_div5(x) x

#define vec_avg3(x0,x1,x2) x0
#define vec_avg5(x0,x1,x2,x3,x4) x0

// subs a with b with abs output
#define vec_subabs(a, b) _mm_subs_epu8(_mm_max_epu8(a,b), _mm_min_epu8(a,b))

// vec compare 
#define vec_cmplsb(a, b, c, d) a = _mm_sub_epi8(a, d);\
							   b = _mm_sub_epi8(b, d);\
							   c = _mm_cmpgt_epi8(a, b);\
							   a = _mm_add_epi8(a, d);\
							   b = _mm_add_epi8(b, d)
#endif // __SIMD_MACRO_H__
