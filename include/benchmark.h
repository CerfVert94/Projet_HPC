#ifndef  __TEST_BENCHMARK_H__
#define __TEST_BENCHMARK_H__

#pragma message("  include  benchmark.h")

struct complete_process_set {
    char func_name[128];
    
    void (*sd_step0)(uint8** X, uint8** Y, uint8** Z, long nrl, long nrh, long ncl, long nch, uint8 n_coeff, uint8 v_min, uint8 v_max);
    void (*sd_func)(p_image t0, p_image t1, uint8 n_coeff, uint8 v_min, uint8 v_max);
    void (*morpho_func)(uint8** X, long nrl, long nrh, long ncl, long nch, uint8 **tempBuffer, uint8 **Y);
    

    void (*vec_sd_step0)(vuint8** X, vuint8** Y, vuint8** Z, long nrl, long nrh, int v0, int v1, uint8 n_coeff, uint8 v_min, uint8 v_max);
    void (*vec_sd_func)(p_vimage t0, p_vimage t1, uint8 n_coeff, uint8 v_min, uint8 v_max);
    void (*vec_morpho_func)(vuint8** vX, int i0, int i1, long ncl, long nch, int j0, int j1, vuint8 **vTempBuffer, vuint8 **vY);
    
    p_image t0;
    p_image t1;
    uint8 **X;
    uint8 **Y;
    uint8 **Z; 
    long nrl;
    long nrh;
    long ncl;
    long nch;
    p_vimage v_t0;
    p_vimage v_t1;
    vuint8 **vX;
    vuint8 **vY;
    vuint8 **vZ; 
    int i0;
    int i1;
    int j0;
    int j1;
    uint8 n_coeff;
    uint8 v_min; 
    uint8 v_max; 
    instruction_type instr_type;
}full_process_set;
void launch_complete_process_benchmark(const char *filename, struct complete_process_set *cps , const int nb_sets, const int nb_tests,const int packet_size, long min_size, long max_size, long step);
void launch_morpho_benchmark          (const char *filename, struct morpho_set *morphos       , const int nb_sets, const int nb_tests,const int packet_size, long min_size, long max_size, long step);
void launch_packed_morpho_benchmark   (const char *filename, struct morpho_set *packed_morphos, const int nb_sets, const int nb_tests,const int packet_size, long min_size, long max_size, long step);
void launch_SD_step_benchmark         (const char *filename, struct sd_set *sd_steps          , const int nb_sets, const int nb_tests,const int packet_size, long min_size, long max_size, long step);
void launch_SD_benchmark              (const char *filename, struct complete_sd_set *csds     , const int nb_sets, const int nb_tests,const int packet_size, long min_size, long max_size, long step);

unsigned long long get_cpu_cycles_of_morpho           (struct morpho_set *ptr_mset                  , uint8 **X,  long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y);
unsigned long long get_min_cpu_cycles_of_morpho       (struct morpho_set *ptr_mset, long packet_size, uint8 **X,  long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y);

unsigned long long get_cpu_cycles_of_packed_morpho    (struct morpho_set *ptr_mset                  , uint8 **X,  long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y);
unsigned long long get_min_cpu_cycles_of_packed_morpho(struct morpho_set *ptr_mset, long packet_size, uint8 **X,  long nrl, long nrh, long ncl, long nch, uint8 **temp_buffer, uint8 **Y);

unsigned long long get_cpu_cycles_of_sd_step          (struct sd_set       *sdset                  , uint8** X, uint8** Y, uint8** Z, long nrl, long nrh, long ncl, long nch, uint8 n_coeff, uint8 v_min, uint8 v_max);
unsigned long long get_min_cpu_cycles_of_sd_step      (struct sd_set       *sdset, long packet_size, uint8** X, uint8** Y, uint8** Z, long nrl, long nrh, long ncl, long nch, uint8 n_coeff, uint8 v_min, uint8 v_max);

unsigned long long get_cpu_cycles_of_sd               (struct complete_sd_set *sdset                  , p_image t0, p_image t1, uint8 n_coeff, uint8 v_min, uint8 v_max);
unsigned long long get_min_cpu_cycles_of_sd           (struct complete_sd_set *sdset, long packet_size, p_image t0, p_image t1, uint8 n_coeff, uint8 v_min, uint8 v_max);


unsigned long long get_cpu_cycles_of_vec_morpho           (struct morpho_set *ptr_mset                  , vuint8 **vX, int i0, int i1, long ncl, long nch, int j0, int j1, vuint8 **temp_vBuffer, vuint8 **vY);
unsigned long long get_min_cpu_cycles_of_vec_morpho       (struct morpho_set *ptr_mset, long packet_size, vuint8 **vX, int i0, int i1, long ncl, long nch, int j0, int j1, vuint8 **temp_vBuffer, vuint8 **vY);

unsigned long long get_cpu_cycles_of_vec_packed_morpho    (struct morpho_set *ptr_mset                  , vuint8 **vX, int i0, int i1, int j0, int j1, vuint8 **temp_vBuffer, vuint8 **vY);
unsigned long long get_min_cpu_cycles_of_vec_packed_morpho(struct morpho_set *ptr_mset, long packet_size, vuint8 **vX, int i0, int i1, int j0, int j1, vuint8 **temp_vBuffer, vuint8 **vY);

unsigned long long get_cpu_cycles_of_vec_sd_step          (struct sd_set       *sdset                  , vuint8 **vX, vuint8 **vY, vuint8 **vZ, int i0, int i1, int j0, int j1, uint8 n_coeff, uint8 v_min, uint8 v_max);
unsigned long long get_min_cpu_cycles_of_vec_sd_step      (struct sd_set       *sdset, long packet_size, vuint8 **vX, vuint8 **vY, vuint8 **vZ, int i0, int i1, int j0, int j1, uint8 n_coeff, uint8 v_min, uint8 v_max);

unsigned long long get_cpu_cycles_of_vec_sd               (struct complete_sd_set *sdset                  , p_vimage t0, p_vimage t1, uint8 n_coeff, uint8 v_min, uint8 v_max);
unsigned long long get_min_cpu_cycles_of_vec_sd           (struct complete_sd_set *sdset, long packet_size, p_vimage t0, p_vimage t1, uint8 n_coeff, uint8 v_min, uint8 v_max);

double **benchmark_of_morpho          (struct morpho_set           *msets , long nb_sets, long ls, long hs, long step, int nb_tests, int packet_size);
double **benchmark_of_packed_morpho   (struct morpho_set           *msets , long nb_sets, long ls, long hs, long step, int nb_tests, int packet_size);
double **benchmark_of_sd_step         (struct sd_set               *sdsets, long nb_sets, long ls, long hs, long step, int nb_tests, int packet_size);
double **benchmark_of_sd              (struct complete_sd_set      *sdsets, long nb_sets, long ls, long hs, long step, int nb_tests, int packet_size);
double **benchmark_of_complete_process(struct complete_process_set *cproc , long nb_sets, long ls, long hs, long step, int nb_tests, int packet_size);

double **init_benchmark_results(long nb_funcs, long size, long step);
void     free_benchmark_results(double **results, long nb_funcs);






#endif /*_BENCHMARK_H_*/