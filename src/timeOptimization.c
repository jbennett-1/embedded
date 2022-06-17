#include "eig_vec_decomp.h"
#include "system.h"
#include "arm_math.h"

#define ARMCM4_FP
#define VEC_NUM 8U //256, 2^8 vectors, bin function after fft to reduce the size (2048->1024)
#define VEC_LEN 256U //length 8

extern float32_t input_data[];

float32_t bin_group[VEC_NUM*8]; //change to 8 with data
float32_t tmp[VEC_LEN/2];
float32_t cov_buffer[(VEC_LEN/2)*(VEC_LEN/2)];
float32_t cov_mat_means[VEC_LEN/2];
float32_t eig_buffer[VEC_LEN/2];
float32_t s_buffer[VEC_LEN];

arm_matrix_instance_f32 matrix_buffer;
arm_rfft_fast_instance_f32 fft_data;

struct fft_pca_args{
    struct eig_decomp_args* eig_args;

    float32_t* eig_buffer;
    float32_t* cov_buffer;
    float32_t* cov_mat_means;
    float32_t* output_buffer;
    float32_t* tmp;
     
    arm_matrix_instance_f32* matrix_buffer;
    uint32_t vec_len;
    uint32_t vec_num;
    float32_t* input_data;	

    arm_rfft_fast_instance_f32* fft_data;
    uint32_t fftLen;
    uint8_t ifftFlag;

};

void initialize(float32_t* s_buffer, float32_t* tmp, struct eig_decomp_args* eig_input, struct fft_pca_args* input, arm_matrix_instance_f32* matrix_buffer, float32_t* input_data, uint32_t vec_len, uint32_t vec_num, arm_rfft_fast_instance_f32* fft_data, uint32_t fftLen, uint8_t ifftFlag, float32_t* cov_buffer, float32_t* cov_mat_means, float32_t* eig_buffer)
{
    input->input_data=input_data;
    input->eig_buffer=eig_buffer;
    input->cov_buffer=cov_buffer;
    input->cov_mat_means=cov_mat_means;

    input->matrix_buffer=matrix_buffer;
    input->vec_len = vec_len;
    input->vec_num = vec_num;

    input->tmp=tmp;
    input->fft_data= fft_data;
    input->fftLen = fftLen;
    input->ifftFlag = ifftFlag;

    eig_input->eig_vec=eig_buffer;
    eig_input->dim_size=vec_len/2;
    eig_input->s=s_buffer;
    eig_input->execs=10000;
    eig_input->err_tol=0.01;
}


void fft_obs_matrix(uint32_t vec_len, uint32_t vec_num, uint8_t ifftFlag, float32_t* input_data)
{
    for (uint32_t i = 0; i < vec_num; i++)
    {
	/*  Since the data is real, the output should be length n/2 + 1. So,
            placing the output as (buffer + i*(vec_len/2 + 1)) allows us to 
            condense the matrix without leaving any extra space because of
            how real ffts work. */
    arm_rfft_fast_f32(&fft_data, (input_data + i*vec_len), tmp, ifftFlag); 
    arm_cmplx_mag_f32(tmp, input_data + i*(vec_len/2), vec_num);
    }
}

float32_t mean(float32_t* vector, uint32_t vec_len, uint32_t vec_num)
{
    float32_t sum = 0;

    for (uint32_t i = 0; i < vec_num; i++)
    {
        sum += vector[i*vec_len]; // because matrices are row major instead of column major
    }

    return sum / vec_num;
}

void calc_means(float32_t* input_data, float32_t* cov_mat_means, uint32_t vec_len, uint32_t vec_num)
{
    for (uint32_t i = 0; i < vec_len; i++)
    {
        cov_mat_means[i] = 0;
        cov_mat_means[i] = mean(&input_data[i], vec_len, vec_num);
    } 
}

void bin(float32_t* input_data, uint32_t vec_len, uint32_t vec_num, uint32_t bins)
{
   uint32_t block_size = vec_len/bins;
   
   for(uint32_t e = 0; e < vec_num; e++) { //iterator through bins
	for(uint32_t i=0; i < bins; i++) { //change to /8 when vec num is changed   
	    for(uint32_t b=0; b < block_size; b++) {
		bin_group[i] = input_data[i*block_size+b]; //change to +7 later
    	    }
	} 
    }
}

void cov(float32_t* input_data, float32_t* cov_buffer, float32_t* cov_mat_means, uint32_t vec_len, uint32_t vec_num)
{
    calc_means(input_data, cov_mat_means, vec_len, vec_num);

    for (int32_t i = 0; i < vec_len; i++) {
        for (int32_t j = i; j < vec_len; j++) {
            // calculate covariance between two vectors
            cov_buffer[i*vec_len + j] = 0;
            for (int32_t k = 0; k < vec_num; k++) {
                cov_buffer[i*vec_len + j] += (input_data[k*vec_len+i]-cov_mat_means[i])*(input_data[k*vec_len+j]-cov_mat_means[j]);
            }
            cov_buffer[i*vec_len + j] /= vec_num - 1;
            if (i != j) {
                cov_buffer[j*vec_len + i] = cov_buffer[i*vec_len + j];
            }
        }
    }
}

void fft_pca(struct fft_pca_args* args, struct eig_decomp_args* eig_input, uint32_t vec_len, uint32_t vec_num, uint8_t ifftFlag)
{ 
    uint32_t bins=8;

    fft_obs_matrix(vec_len, vec_num, ifftFlag,args->input_data);   
    vec_len = (vec_len/2); // since data is real, vectors after fft are length n/2 + 1     input_data = data_buffer
    bin(args->input_data,vec_len,vec_num, bins);

    cov(args->input_data, cov_buffer, cov_mat_means, vec_len, vec_num);
    
    eig_decomp(&matrix_buffer, eig_input);
}

void main(){
    struct fft_pca_args args; 
    struct eig_decomp_args eig_input;
    extern uint32_t tick;
    tick=0;

    uint8_t ifftFlag=0;
    uint32_t vec_len = VEC_LEN;
    uint32_t vec_num = VEC_NUM;
    uint32_t fftLen=VEC_LEN;

    initialize(s_buffer, tmp, &eig_input, &args, &matrix_buffer, input_data, vec_len, vec_num, &fft_data, fftLen, ifftFlag, cov_buffer, cov_mat_means, eig_buffer);

    SysTick->CTRL=0;
    SysTick->VAL=0;
    SysTick_Config((SystemCoreClock/500)); //parameter becomes reload value
    NVIC_EnableIRQ(SysTick_IRQn); 
    asm("cpsie if");

    arm_rfft_fast_init_f32(&fft_data, fftLen);  
    arm_mat_init_f32(&matrix_buffer, vec_len, vec_num, input_data);

    fft_pca(&args, &eig_input, vec_num, vec_len, ifftFlag);
}

