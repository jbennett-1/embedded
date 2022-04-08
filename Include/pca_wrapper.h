#ifndef _PCA_WRAPPER_H
#define _PCA_WRAPPER_H

#include "headers.h"

struct fft_pca_args {
    struct fft_args* f_args;
    uint32_t vec_len;
    uint32_t vec_num;
    float* fft_buf;
    float* cov_mat;
    float* cov_mat_means;
    float* princ_comp;
    // Extra memory for ssyevx
    uint32_t* iwork;
    float* work;
    uint32_t lwork;
};

void fft_pca(float* input_buffer, float* output_buffer, void* args);
struct fft_pca_args* alloc_fft_pca_args(uint32_t vec_len, uint32_t vec_num);
void free_fft_pca_args(struct fft_pca_args* args);

#endif
