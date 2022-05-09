#include <stdint.h>
#include "eig_vec_decomp.h"
#include "arm_math.h"
#define ARMCM7

struct fft_pca_args{
	arm_matrix_instance_f32* input_buffer;
	float32_t* eig_buffer;
	float32_t* cov_buffer;
	float32_t* output_buffer;

	arm_rfft_fast_instance_f32* fft_data;
	
	float32_t* cov_mat;
	float32_t* cov_mat_means;
	
	uint16_t vec_len;
	uint16_t vec_num;
	
	uint8_t ifftFlag;
	uint8_t bits;
};

void initialize(void* args, void *eig_args, float32_t* cov_buffer, float32_t* eig_buffer, float32_t* output_buffer, float32_t* fft_input, uint8_t ifftFlag, uint8_t bits, uint16_t vec_len, uint16_t vec_num, float32_t* cov_mat_means)
{
    struct fft_pca_args* input = (struct fft_pca_args*) args;
//    struct eig_decomp_args* eig_input = (struct eig_decomp_args*) eig_args;
   
    //do I need to initialize this?
    //eig_input->eig_args=eig_args;
    
    input->cov_buffer=cov_buffer;
    input->eig_buffer=eig_buffer;
    input->fft_data=(arm_rfft_fast_instance_f32*) fft_input;

    input->vec_len = vec_len;
    input->vec_num = vec_num;

    input->ifftFlag = ifftFlag;
    input->bits = bits;

    input->cov_mat=cov_buffer;

    input->cov_mat_means=cov_mat_means;
   
    arm_rfft_fast_init_f32(input->fft_data, vec_len);
   
    arm_mat_init_f32(input->input_buffer, vec_len, vec_num, data);
}

void fft_obs_matrix(struct fft_pca_args* args)
{
    int16_t vec_num=args->vec_num;
    int16_t vec_len=args->vec_len;
	
    float32_t *p = args->data;
    float32_t *pOut = (float32_t*) args->input_buffer;

    for (int16_t i = 0; i < vec_num; i++)
    {
	/*  Since the data is real, the output should be length n/2 + 1. So,
            placing the output as (buffer + i*(vec_len/2 + 1)) allows us to 
            condense the matrix without leaving any extra space because of
            how real ffts work. */
    arm_rfft_fast_f32(args->fft_data,p + i*vec_len, pOut+ i*(vec_num/2+1), args->ifftFlag);
    }	
}

void cov(struct fft_pca_args* args)
//float32_t* A, float32_t* cov_mat, float32_t* means, uint32_t vec_len, uint32_t vec_num)
{
//in pca code, float A = input_buffer, which corresponds to fft_buff in this code
//will this syntax not access the struct member but instead take the value from it as a copy
    float32_t* inp_buf = args->data;
    uint32_t vec_len = args->vec_len;
    uint32_t vec_num = args->vec_num;
    float32_t* means = args->cov_mat_means;
    float32_t* cov_mat = args->cov_mat;

    for (int16_t i = 0; i < vec_len; i++) {
        for (int16_t j = i; j < vec_len; j++) {
            // calculate covariance between two vectors
            cov_mat[i*vec_len + j] = 0;
            for (int16_t k = 0; k < vec_num; k++) {
                cov_mat[i*vec_len + j] += (inp_buf[k*vec_len + i] - means[i]) * (inp_buf[k*vec_len + j] - means[j]);
            }
            cov_mat[i*vec_len + j] /= vec_num - 1;
            if (i != j) {
                cov_mat[j*vec_len + i] = cov_mat[i*vec_len + j];
            }
        }
    }
}

void fft_pca(struct fft_pca_args* args)
{
   uint16_t vec_len = args->vec_len;
   struct eig_decomp_args* decomp_args;

   fft_obs_matrix(args);

   vec_len = (vec_len/2 + 1); // since data is real, vectors after fft are length n/2 + 1

   cov(args);
   eig_decomp(decomp_args);
}

void main(){
   float32_t input[] = {1, 2, 6, 4, 3, 9, 7, 0, 12, 4, 7, 23, 9, 2, 7 , 3, 1, 3, 7, 1, 2, 9, 2, 1, 4, 5, 3, 1, 8, 2, 1, 1};
   float32_t output[] = {0, 0};

   struct fft_pca_args* args;
   struct eig_decomp_args* decomp_args;

   uint8_t ifftFlag = 0;
   uint8_t bits = 4;
   uint16_t vec_len = 16;
   uint16_t vec_num = 2;
   
   float32_t cov_buff[vec_len*vec_len];
   float32_t* cov_mat_means=(float32_t) vec_len;
   float32_t eig_buff[32];

   initialize((void*) args,(void*) decomp_args,cov_buff, eig_buff, output, input, ifftFlag, bits, vec_len, vec_num, cov_mat_means);
}

