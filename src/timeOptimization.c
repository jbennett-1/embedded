#include "pca_wrapper.h"
#include "f2c.h"
#include "clapack.h"
#include "fft_wrapper.h"

void fft_pca(float* input_buffer, float* output_buffer, void* args)
{
    struct fft_pca_args* inputs = (struct fft_pca_args*) args;
    float* fft_buf = inputs->fft_buf;
    float* cov_mat = inputs->cov_mat;
    float* cov_mat_means = inputs->cov_mat_means;
    float* princ_comp = inputs->princ_comp;
    float* work = inputs->work;
    uint32_t* iwork = inputs->iwork;
    uint32_t lwork = inputs->lwork;
    uint32_t vec_len = inputs->vec_len;
    uint32_t vec_num = inputs->vec_num;

    fft_obs_matrix(input_buffer, fft_buf,  vec_len, vec_num, inputs->f_args);

    vec_len = (vec_len/2 + 1); // since data is real, vectors after fft are length n/2 + 1 

    cov(fft_buf, cov_mat, cov_mat_means, vec_len, vec_num);

    integer n = vec_len, lda = vec_len, ldz = vec_len;
    integer lower_limit = vec_len;
    integer upper_limit = vec_len;
    integer eig_val_found = 0;
    integer ifail = 0;
    integer info = 0;
    real eig_val = 0;
    real absol = 0.001;
//   ssyevx_("Vectors", "Indices", "Upper", (integer*) &n, (real*) cov_mat, (integer*) &lda, (real*) 0, (real*) 0, (integer*) &upper_limit, (integer*) &lower_limit,

//(real*) &absol, (integer*) &eig_val_found, (real*) &eig_val, (real*) princ_comp, (integer*) &ldz, (real*) work, (integer*) &lwork, (integer*) iwork,
  //         (integer*) &ifail, (integer*) &info);

   // fix_output(princ_comp, vec_len);

   // project_data(fft_buf, princ_comp, output_buffer, vec_len, vec_num);
}


void fft_obs_matrix(float* input, float* output, uint32_t vec_len, uint32_t vec_num, struct fft_args* args)
{
    for (int i = 0; i < vec_num; i++)
    {
        /*
            Since the data is real, the output should be length n/2 + 1. So,
            placing the output as (buffer + i*(vec_len/2 + 1)) allows us to 
            condense the matrix without leaving any extra space because of
            how real ffts work.
        */
        fft_abs((input + i*vec_len), (output + i*(vec_len/2 + 1)), args);
    }
}


void cov(float* A, float* cov_mat, float* means, uint32_t vec_len, uint32_t vec_num)
{

    calc_means(A, means, vec_len, vec_num);

    for (int i = 0; i < vec_len; i++) {

        for (int j = i; j < vec_len; j++) {
            // calculate covariance between two vectors
            cov_mat[i*vec_len + j] = 0;
            for (int k = 0; k < vec_num; k++) {
                cov_mat[i*vec_len + j] += (A[k*vec_len + i] - means[i]) * (A[k*vec_len + j] - means[j]);
            }
            cov_mat[i*vec_len + j] /= vec_num - 1;
            if (i != j) {
                cov_mat[j*vec_len + i] = cov_mat[i*vec_len + j];
            }
        }
    }
}

