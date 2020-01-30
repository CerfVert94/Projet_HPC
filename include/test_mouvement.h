/* ------------------------ */
/* --- test_mouvement.h --- */
/* ------------------------ */
struct sd_set{
    char func_name[128];
    void (*sd_func)(uint8**, uint8**, uint8**, long, long, long, long);
    enum {NO_PACK, HPACK, VPACK}pack_type;
};


uint8 test_corps_SigmaDelta_step1(uint8 t_1M, uint8 tI);
uint8 test_corps_SigmaDelta_step2(uint8 tM, uint8 tI);
uint8 test_corps_SigmaDelta_step3(uint8 t_1V, uint8 tO);
uint8 test_corps_SigmaDelta_step4(uint8 tO, uint8 tV);

void all_test();