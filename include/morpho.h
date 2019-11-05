typedef struct struct_elem {
	long orix;
	long oriy;
	long nrow;
	long ncol;
	uint8** m;
} struct_elem;

typedef struct struct_elem struct_elem, *p_struct_elem;