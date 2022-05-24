#include <stdint.h>
#include "eig_vec_decomp.h"
#include "system.h"
#include "arm_math.h"
#define ARMCM4_FP
#define VEC_NUM 4U
#define VEC_LEN 8U
#define FFT_LEN 32U

float32_t input_data[VEC_LEN*VEC_NUM];
float32_t output_buffer[VEC_LEN*VEC_NUM];
float32_t tmp[VEC_LEN/2+1];
float32_t pTwiddleRFFT[FFT_LEN];

float32_t cov_buffer[(VEC_LEN/2+1)*(VEC_LEN/2+1)];
float32_t cov_mat_means[VEC_LEN/2+1];
float32_t eig_buffer[VEC_LEN/2+1];

struct fft_pca_args{
    struct eig_decomp_args* eig_decomp;

    float32_t* eig_buffer;
    float32_t* cov_buffer;
    float32_t* cov_mat_means;
    float32_t* output_buffer;

    arm_matrix_instance_f32* matrix_buffer;
    uint16_t vec_len;
    uint16_t vec_num;
    float32_t* input_data;	

    arm_rfft_fast_instance_f32* fft_data;
    float32_t* pTwiddleRFFT;
    uint16_t fftLen;
    arm_cfft_instance_f32* cFFT;
    uint8_t ifftFlag;

    uint8_t bits;
};

void initialize(float32_t* pTwiddleRFFT, arm_cfft_instance_f32* cFFT, float32_t* tmp, struct eig_decomp_args* eig_input, struct fft_pca_args* input, float32_t* output_buffer, arm_matrix_instance_f32* matrix_buffer, float32_t* input_data, uint16_t vec_len, uint16_t vec_num, arm_rfft_fast_instance_f32* fft_data, uint16_t fftLen, uint8_t ifftFlag, float32_t* cov_buffer, float32_t* cov_mat_means, float32_t* eig_buffer, uint8_t bits)
{
    input->eig_buffer=eig_buffer;
    input->cov_buffer=cov_buffer;
    input->cov_mat_means=cov_mat_means;
    input->output_buffer=output_buffer;
    input->input_data=input_data;

    input->cFFT=cFFT;
    input->pTwiddleRFFT=pTwiddleRFFT;

    input->matrix_buffer=matrix_buffer;
    input->vec_len = vec_len;
    input->vec_num = vec_num;

    input->fft_data= fft_data;
    input->fftLen = fftLen;

    input->ifftFlag = ifftFlag;
    input->bits = bits;

    eig_input->eig_vec=eig_buffer;

    arm_rfft_fast_init_f32(fft_data, fftLen);
    arm_mat_init_f32(matrix_buffer, vec_len, vec_num, input_data);
}

void fft_obs_matrix(uint16_t fftLen, uint8_t ifftFlag, float32_t* input_data, arm_rfft_fast_instance_f32* fft_data, float32_t* tmp, uint16_t vec_len, uint16_t vec_num, arm_cfft_instance_f32* cFFT)
{
    arm_status status;
    status=ARM_MATH_SUCCESS;
    status=arm_cfft_init_f32(cFFT,fftLen);

    for (int16_t i = 0; i < vec_num; i++)
    {
	/*  Since the data is real, the output should be length n/2 + 1. So,
            placing the output as (buffer + i*(vec_len/2 + 1)) allows us to 
            condense the matrix without leaving any extra space because of
            how real ffts work. */
    arm_rfft_fast_f32(fft_data, (input_data + i*vec_len), tmp, ifftFlag);
    arm_cmplx_mag_f32(tmp, input_data + i*(vec_len/2+1), vec_num);
    }
}

void cov(float32_t* input_data, float32_t* cov_buffer, float32_t* cov_mat_means, uint16_t vec_len, uint16_t vec_num)
{
    for (int16_t i = 0; i < vec_len; i++) {
        for (int16_t j = i; j < vec_len; j++) {
            // calculate covariance between two vectors
            cov_buffer[i*vec_len + j] = 0;
            for (int16_t k = 0; k < vec_num; k++) {
                cov_buffer[i*vec_len + j] += (input_data[k*vec_len+i]-cov_mat_means[i])*(input_data[k*vec_len+j]-cov_mat_means[j]);
            }
            cov_buffer[i*vec_len + j] /= vec_num - 1;
            if (i != j) {
                cov_buffer[j*vec_len + i] = cov_buffer[i*vec_len + j];
            }
        }
    }
}

void fft_pca(struct fft_pca_args* args)
{
    uint16_t vec_len=args->vec_len; 
    uint16_t vec_num=args->vec_num; 
    uint16_t fftLen=args->fftLen;
    uint8_t ifftFlag = args->ifftFlag;

    arm_rfft_fast_instance_f32* fft_data = args->fft_data;
    arm_matrix_instance_f32* data_matrix = args->matrix_buffer;
    arm_cfft_instance_f32* cFFT=args->cFFT;

    fft_obs_matrix(ifftFlag, fftLen, input_data, fft_data, tmp, vec_len, vec_num, cFFT);
   
    vec_len = (vec_len/2 + 1); // since data is real, vectors after fft are length n/2 + 1     input_data = data_buffer
    
    cov(input_data, args->cov_buffer, args->cov_mat_means, vec_len, vec_num);
    
    eig_decomp(data_matrix, args->eig_decomp);
}

void main(){
    struct fft_pca_args args; 
    struct eig_decomp_args eig_input;
   
    arm_matrix_instance_f32 matrix_buffer;
    arm_rfft_fast_instance_f32 fft_data;
    arm_cfft_instance_f32 cFFT;    

    uint16_t vec_len = VEC_LEN;
    uint16_t vec_num = VEC_NUM;
    uint16_t fftLen = FFT_LEN;

    uint8_t ifftFlag = 0;
    uint8_t bits = 4;

    initialize(pTwiddleRFFT, &cFFT, tmp, &eig_input, &args, output_buffer, &matrix_buffer, input_data, vec_len, vec_num, &fft_data, fftLen, ifftFlag, cov_buffer, cov_mat_means, eig_buffer, bits);
    fft_pca(&args);
}

