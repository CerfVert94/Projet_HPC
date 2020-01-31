/* ------------------------ */
/* --- test_mouvement.h --- */
/* ------------------------ */

uint8 test_corps_SigmaDelta_step1(uint8 t_1M, uint8 tI);
uint8 test_corps_SigmaDelta_step2(uint8 tM, uint8 tI);
uint8 test_corps_SigmaDelta_step3(uint8 t_1V, uint8 tO);
uint8 test_corps_SigmaDelta_step4(uint8 tO, uint8 tV);

void all_test();

void test_implementation_SD_complete();
void test_implementation_SigmaDelta_step1(struct sd_set *sd, int nb_imps);
void test_implementation_SigmaDelta_step2(struct sd_set *sd, int nb_imps);
void test_implementation_SigmaDelta_step3(struct sd_set *sd, int nb_imps);
void test_implementation_SigmaDelta_step4(struct sd_set *sd, int nb_imps);
bool SD_step1_produces_valid_output(uint8 m_t0,          uint8 i_t1,  uint8 m_t1);
bool SD_step2_produces_valid_output(uint8 m_t ,          uint8 i_t ,  uint8 o_t);
bool SD_step3_produces_valid_output(uint8 v_t0, uint8 n, uint8 o_t, uint8 v_t1, uint8 _vmin, uint8 _vmax);
bool SD_step4_produces_valid_output(uint8 o_t ,          uint8 v_t, uint8 e);