#include <stdint.h>
#include "eig_vec_decomp.h"
#include "arm_math.h"
#define ARMCM7_DP
#define VEC_NUM 2
#define VEC_LEN 16

float32_t input_data[VEC_LEN*VEC_NUM];
float32_t output_buffer[VEC_LEN*VEC_NUM];
float32_t tmp[VEC_LEN/2+1];
float32_t data[VEC_LEN*VEC_NUM];

float32_t cov_buffer[(VEC_LEN/2+1)*(VEC_LEN/2+1)];
float32_t cov_mat_means[VEC_LEN/2+1];
float32_t eig_buffer[VEC_LEN/2+1];
float32_t fft_out[VEC_LEN*VEC_NUM];

struct fft_pca_args{
    struct eig_decomp_args* eig_decomp;

    float32_t* eig_buffer;
    float32_t* cov_buffer;
    float32_t* cov_mat_means;
    float32_t* output_buffer;

    arm_matrix_instance_f32* matrix_buffer;
    uint16_t vec_len;
    uint16_t vec_num;
    float32_t* data;	

    arm_rfft_fast_instance_f32* fft_data;

    uint16_t fft_len;
    float32_t* fft_out;
    uint8_t ifftFlag;

    uint8_t bits;
};

void initialize(float32_t* tmp, float32_t* data, struct eig_decomp_args* eig_input, struct fft_pca_args* input, float32_t* output_buffer, arm_matrix_instance_f32* matrix_buffer, float32_t* input_data, uint16_t vec_len, uint16_t vec_num, arm_rfft_fast_instance_f32* fft_data, uint16_t fft_len, float32_t* fft_out, uint8_t ifftFlag, float32_t* cov_buffer, float32_t* cov_mat_means, float32_t* eig_buffer, uint8_t bits)
{
    input->eig_buffer=eig_buffer;
    input->cov_buffer=cov_buffer;
    input->cov_mat_means=cov_mat_means;
    input->output_buffer=output_buffer;
    input->data=data;

    input->matrix_buffer=matrix_buffer;
    input->vec_len = vec_len;
    input->vec_num = vec_num;

    input->fft_data= fft_data;
    input->fft_len = fft_len;
    input->fft_out = fft_out;

    input->ifftFlag = ifftFlag;
    input->bits = bits;

    eig_input->eig_vec=eig_buffer;

    arm_rfft_fast_init_f32(fft_data, fft_len);
    arm_mat_init_f32(matrix_buffer, vec_len, vec_num, input_data);
}

void fft_obs_matrix(float32_t* data, arm_rfft_fast_instance_f32* fft_data, uint16_t vec_len, uint16_t vec_num, struct fft_pca_args* args)
{
    for (int16_t i = 0; i < vec_num; i++)
    {
	/*  Since the data is real, the output should be length n/2 + 1. So,
            placing the output as (buffer + i*(vec_len/2 + 1)) allows us to 
            condense the matrix without leaving any extra space because of
            how real ffts work. */
    arm_rfft_fast_f32(fft_data, data + i*vec_len, args->fft_out+ i*(vec_num/2+1), args->ifftFlag);
    }	
}

void cov(float32_t* input_data, float32_t* cov_mat, float32_t* means, uint16_t vec_len, uint16_t vec_num)
{
    for (int16_t i = 0; i < vec_len; i++) {
        for (int16_t j = i; j < vec_len; j++) {
            // calculate covariance between two vectors
            cov_mat[i*vec_len + j] = 0;
            for (int16_t k = 0; k < vec_num; k++) {
                cov_mat[i*vec_len + j] += (input_data[k*vec_len + i] - means[i]) * (input_data[k*vec_len + j] - means[j]);
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
    uint16_t vec_len=args->vec_len; 
    uint16_t vec_num=args->vec_num; 
    
    arm_matrix_instance_f32 mat = args->eig_input->targ_mat;
    arm_rfft_fast_instance_f32 fft_data = args->fft_data;
    float32_t input_data = args->data;

    fft_obs_matrix(&input_data, &fft_data, vec_len, vec_num, args);
   
    vec_len = (vec_len/2 + 1); // since data is real, vectors after fft are length n/2 + 1
    
    cov(&input_data, args->cov_buffer, args->cov_mat_means, vec_len, vec_num);
   
    float32_t* eig_input = args->eig_decomp->eig_vec;  
    
    eig_decomp(&mat, eig_input);
    //fix_output(eig_vec, vec_len)
    //project_data(fft_buf, eig_vec, output_buffer, vec_len, vec_num);
}

void main(){
    struct fft_pca_args args; 
    struct eig_decomp_args eig_input;
   
    arm_matrix_instance_f32 matrix_buffer;
    arm_rfft_fast_instance_f32 fft_data;
    
    uint16_t vec_len = VEC_LEN;
    uint16_t vec_num = VEC_NUM;
    uint16_t fft_len = VEC_LEN;

    uint8_t ifftFlag = 0;
    uint8_t bits = 4;

    initialize(tmp, data, &eig_input, &args, output_buffer, &matrix_buffer, input_data, vec_len, vec_num, &fft_data, fft_len, fft_out, ifftFlag, cov_buffer, cov_mat_means, eig_buffer, bits);
    fft_pca(&args);
}

