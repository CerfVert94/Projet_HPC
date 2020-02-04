/* -------------------- */
/* --- simd_macro.h --- */
/* -------------------- */


#ifndef __SIMD_MACRO_H__
#define __SIMD_MACRO_H__

#pragma message("  include  simd_macro.h")

// -------
// defines
// -------

#define leftshift1a  _mm_setr_epi8(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,-1)
#define leftshift1b  _mm_setr_epi8(-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0)

#define rightshift1a _mm_setr_epi8(-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14)
#define rightshift1b _mm_setr_epi8(15,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1)


#define leftshift2a  _mm_setr_epi8(2,3,4,5,6,7,8,9,10,11,12,13,14,15,-1,-1)
#define leftshift2b  _mm_setr_epi8(-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,1)

#define rightshift2a _mm_setr_epi8(-1,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13)
#define rightshift2b _mm_setr_epi8(14,15,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1)




// -----------------------
// melanges / permutations
// -----------------------

#define vec_border(a,b) b
#define vec_middle(a,b) a

#define vec_left1(v0, v1) (_mm_or_si128(_mm_shuffle_epi8(v0, leftshift1a), _mm_shuffle_epi8(v1, leftshift1b)))
#define vec_left2(v0, v1) (_mm_or_si128(_mm_shuffle_epi8(v0, leftshift2a), _mm_shuffle_epi8(v1, leftshift2b)))
#define vec_left3(v0, v1) v1
#define vec_left4(v0, v1) v1

#define vec_right1(v1, v2) (_mm_or_si128(_mm_shuffle_epi8(v2, rightshift1a), _mm_shuffle_epi8(v1, rightshift1b)))
#define vec_right2(v1, v2) (_mm_or_si128(_mm_shuffle_epi8(v2, rightshift2a), _mm_shuffle_epi8(v1, rightshift2b)))
#define vec_right3(v1, v2) v1
#define vec_right4(v1, v2) v1 

// ----------
// operateurs
// ----------

#define vector_and3(A, B, C)	       (_mm_and_si128(_mm_and_si128(A, B), C))
#define vector_and3_row1shift(A, B, C) (vector_and3(vec_right1(A, B), B, vec_left1(B, C)))
#define vector_and3_row2shift(A, B, C) (vector_and3(vec_right2(A, B), B, vec_left2(B, C)))
#define vector_and5(A, B, C, D, E)     (vector_and3(vector_and3(A,B,C), D, E))

#define vector_or3(A, B, C)	    	   (_mm_or_si128(_mm_or_si128(A, B),C))
#define vector_or3_row1shift(A, B, C)  (vector_or3(vec_right1(A, B), B, vec_left1(B, C)))
#define vector_or3_row2shift(A, B, C)  (vector_or3(vec_right2(A, B), B, vec_left2(B, C)))
#define vector_or5(A, B, C, D, E)      (vector_or3(vector_or3(A,B,C), D, E))

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
#define vec_cmpgt(a, b, c, d) a = _mm_sub_epi8(a, d);\
							   b = _mm_sub_epi8(b, d);\
							   c = _mm_cmpgt_epi8(a, b);\
							   a = _mm_add_epi8(a, d);\
							   b = _mm_add_epi8(b, d)

#define vec_cmplt(a, b, c, d) a = _mm_sub_epi8(a, d);\
							   b = _mm_sub_epi8(b, d);\
							   c = _mm_cmplt_epi8(a, b);\
							   a = _mm_add_epi8(a, d);\
							   b = _mm_add_epi8(b, d)

#define vec16_cmplt(a, b, c, d) a = _mm_sub_epi16(a, d);\
							   b = _mm_sub_epi16(b, d);\
							   c = _mm_cmplt_epi16(a, b);\
							   a = _mm_add_epi16(a, d);\
							   b = _mm_add_epi16(b, d)


#endif // __SIMD_MACRO_H__
