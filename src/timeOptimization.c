#include <stdint.h>
#include "dsp/transform_functions.h"

//for headers, replace anything fft with the libraries
void fft_pca(float* input_buffer, float* output_buffer, void* args)
{
	//keep args in functions, args
	int* input_buf[] = {1, 2, 6, 4, 3, 9, 7, 0, 12, 4, 7, 23, 9, 2, 7 , 3, 1, 3, 7, 1, 2, 9, 2, 1, 4, 5, 3, 1, 8, 2, 1, 1};
	struct fft_pca_args* inputs = (struct fft_pca_args*) args;
	
	//data
//	float* fft_buf[] = *output_buffer;
	//output buffer
	float* fft_buf = inputs->fft_buf;
	float* cov_mat = inputs->cov_mat;
	float* cov_mat_means = inputs->cov_mat_means;
	uint32_t vec_len = inputs->vec_len;
	uint32_t vec_num = inputs->vec_num;
	
	fft_obs_matrix(input_buf, fft_buf,  vec_len, vec_num);
	vec_len = (vec_len/2 + 1); // since data is real, vectors after fft are length n/2 + 1 
	cov(fft_buf, cov_mat, cov_mat_means, vec_len, vec_num);
}


void fft_obs_matrix(float* input, float* output, uint32_t vec_len, uint32_t vec_num,struct fft_args* args)
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
}
/*
//not specified for code, rework this section out
struct fft_pca_args* alloc_fft_pca_args(uint32_t vec_len, uint32_t vec_num)
{
    struct fft_pca_args* args = malloc(sizeof(struct fft_pca_args));
    
    args->f_args = alloc_fft_args(vec_len);
    args->vec_len = vec_len;
    args->vec_num = vec_num;
    args->fft_buf = malloc(sizeof(float)*(vec_len/2 + 1)*vec_num);
    args->cov_mat = malloc(sizeof(float)*vec_len*vec_len);
    args->cov_mat_means = malloc(sizeof(float)*vec_len);

    // Allocate correct amount of memory for work
    int lower_limit = vec_len;
    int upper_limit = vec_len;
    args->work = malloc(sizeof(float)*50000);

    return args;
}
*/
//overun buffer before going through fft, adjust for when the buffer changes after it halves
