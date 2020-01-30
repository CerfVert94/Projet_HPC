#ifndef  __TEST_MORPHO_H__
#define __TEST_MORPHO_H__

#pragma message("  include  test_morpho.h")

// A struct type that contains a morpho function and a structuring element
struct morpho_set{
    char func_name[128];
    void (*morpho_func)(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y);
    enum {NO_PACK, HPACK, VPACK}pack_type;
};

struct sd_set{
    char func_name[128];
    void (*sd_func)(uint8** X, uint8** Y, uint8** Z, long nrl, long nrh, long ncl, long nch);
};



void test_dilations (struct morpho_set *dilation_sets, const int nb_implementations, bool display);
void test_erosions  (struct morpho_set *erosion_sets , const int nb_implementations, bool display);
void test_sequences (struct morpho_set *sequence_sets, const int nb_implementations, bool display);

void test_implementation_erosion3  (struct morpho_set *morpho_set);
void test_implementation_dilation3 (struct morpho_set *morpho_set);
void test_implementation_erosion5  (struct morpho_set *morpho_set );
void test_implementation_dilation5 (struct morpho_set *morpho_set );

void test_intergration(char *filename, struct morpho_set *naive_morpho_set, struct morpho_set *morpho_sets, const int nb_implementations, bool display);
void test_packed_intergration(char *filename, struct morpho_set *naive_morpho_set, struct morpho_set *morpho_sets, const int nb_implementations, bool display);

// void test_integration_erosion_5x5 (struct morpho_set *erosion_set);
// void test_integration_dilation_5x5(struct morpho_set *dilation_set);

bool morpho_produces_one(struct morpho_set *mset, uint8** X);

bool check_dimension_of_square_structuring_element(p_struct_elem_dim s,  long size);


static inline uint32 extract_bits_from_permutation(uint32 perm, long irow, long ncol);
static inline uint32 extract_bits_from_permutation(uint32 perm, long irow, long ncol) {
    return (perm >> irow * ncol);
}

static inline uint8 get_column_at(uint8 col_vals, long icol);
static inline uint8 get_column_at(uint8 col_vals, long icol) {
    return (col_vals >> icol) & 0x1;
}

uint8 **rand_ui8matrix(long size, p_struct_elem_dim s);
uint8 **ui8matrix_checker(long nrl, long nrh, long ncl, long nch, const long chkr_size, const uint8 xor_mask);
uint8**  ui8matrix_permutation (uint8** m, long nrl, long nrh, long ncl, long nch, uint32 perm);


#endif /*  __TEST_MORPHO_H__ */