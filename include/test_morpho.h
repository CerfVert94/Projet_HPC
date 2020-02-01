#ifndef  __TEST_MORPHO_H__
#define __TEST_MORPHO_H__

#pragma message("  include  test_morpho.h")



void test_dilations(char *filename, struct morpho_set *dilation_sets, const int nb_sets, bool logging);
void test_erosions (char *filename, struct morpho_set *erosion_sets , const int nb_sets, bool logging);
void test_sequences(char *filename, struct morpho_set *sequence_sets, const int nb_sets, bool logging);

void test_implementation_erosion3  (struct morpho_set *morpho_set);
void test_implementation_dilation3 (struct morpho_set *morpho_set);
void test_implementation_erosion5  (struct morpho_set *morpho_set);
void test_implementation_dilation5 (struct morpho_set *morpho_set);

void test_intergration(uint8 **image, long nrl, long nrh, long ncl, long nch, const char *filename, struct morpho_set *naive_morpho_set, struct morpho_set *morpho_sets, int nb_sets, bool logging);
void test_packed_intergration(char *filename, struct morpho_set *naive_morpho_set, struct morpho_set *morpho_sets, const int nb_implementations, bool logging);

// void test_integration_erosion_5x5 (struct morpho_set *erosion_set);
// void test_integration_dilation_5x5(struct morpho_set *dilation_set);

bool morpho_produces_one(struct morpho_set *mset, uint8** X);




static inline uint32 extract_bits_from_permutation(uint32 perm, long irow, long ncol);
static inline uint32 extract_bits_from_permutation(uint32 perm, long irow, long ncol) {
    return (perm >> irow * ncol);
}

static inline uint8 get_column_at(uint8 col_vals, long icol);
static inline uint8 get_column_at(uint8 col_vals, long icol) {
    return (col_vals >> icol) & 0x1;
}

// uint8 **rand_ui8matrix(long size);
uint8 **ui8matrix_checker(long nrl, long nrh, long ncl, long nch, const long chkr_size, const uint8 xor_mask);
uint8**  ui8matrix_permutation (uint8** m, long nrl, long nrh, long ncl, long nch, uint32 perm);


#endif /*  __TEST_MORPHO_H__ */