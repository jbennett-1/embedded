#include <stdint.h>
#include "eig_vec_decomp.h"
#include "system.h"
#include "arm_math.h"
#include "arm_const_structs.h"

#define ARMCM4_FP
#define VEC_NUM 2U
#define VEC_LEN 256U


float32_t input_data[]={0,0,0,0, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256,257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 415, 416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450, 451, 452, 453, 454, 455, 456, 457, 458, 459, 460, 461, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 475, 476, 477, 478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499, 500, 501, 502, 503, 504, 505, 506, 507, 508, 0,0,0,0};
float32_t tmp[VEC_LEN/2];
float32_t cov_buffer[(VEC_LEN/2)*(VEC_LEN/2)];
float32_t cov_mat_means[VEC_LEN/2];
float32_t eig_buffer[VEC_LEN/2];
float32_t s_buffer[VEC_LEN];

static volatile unsigned int tick=0;

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
    eig_input->execs=100;
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
    fft_obs_matrix(vec_len, vec_num, ifftFlag,args->input_data);
   
    vec_len = (vec_len/2); // since data is real, vectors after fft are length n/2 + 1     input_data = data_buffer
    
    cov(args->input_data, cov_buffer, cov_mat_means, vec_len, vec_num);
    
    eig_decomp(&matrix_buffer, eig_input);
}
void SysTick_Handler(){
    tick++;
    if(SysTick->LOAD==0)
	return tick;
}

void main(){
    struct fft_pca_args args; 
    struct eig_decomp_args eig_input;
    
    uint8_t ifftFlag=0;
    uint32_t vec_len = VEC_LEN;
    uint32_t vec_num = VEC_NUM;
    uint32_t fftLen=VEC_LEN;
    uint32_t start_time, stop_time, cycle_count;

    SysTick->CTRL=0;

    SysTick_Config(SystemCoreClock / 1000);
    start_time=(SystemCoreClock/1000) - (SysTick->VAL);

    initialize(s_buffer, tmp, &eig_input, &args, &matrix_buffer, input_data, vec_len, vec_num, &fft_data, fftLen, ifftFlag, cov_buffer, cov_mat_means, eig_buffer);

    arm_rfft_fast_init_f32(&fft_data, fftLen);
    arm_mat_init_f32(&matrix_buffer, vec_len, vec_num, input_data);

    fft_pca(&args, &eig_input, vec_len, vec_num, ifftFlag);
   // stop_time=SysTick->VAL;
    //cycle_count = start_time - stop_time;
}

