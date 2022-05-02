#include <stdint.h>
#include "eig_vec_decomp.h"
#include "arm_math.h"
#define ARMCM7

struct fft_pca_args{
	struct eig_decomp_args *eig_args;
	arm_matrix_instance_f32* input_buffer;
	float32_t* eig_buffer;
	float32_t* cov_buffer;
	float32_t* output_buffer;

	arm_rfft_fast_instance_f32 fft_buffer;
	float32_t* data;
	float32_t* cov_mat_means;
	
	uint16_t vec_len;
	uint16_t vec_num;
	
	uint8_t ifftFlag;
	uint8_t bits;
};

//dynamic mem management//calloc() but write your own version//alloc at compile time//constants in ld for mem 
//need 4 buffers, data, cov, eigen vector buffer + store the product in output buffer
/* matrix is stored as an array, contigious memory in stack, need pointer array (for the buffers of elements x1-x2)
//using headerfiles definitions of arm_matrix instance etc etc
*/

void initialize(struct fft_pca_args* args, struct eig_decomp_args *eig_args, float32_t* cov_buffer, float32_t* eig_buffer, float32_t* output_buffer, arm_rfft_fast_instance_f32* fft_buffer, uint8_t ifftFlag, uint8_t bits, uint16_t vec_len, uint16_t vec_num, float32_t* data, float32_t* cov_mat_means)
{

   args->eig_args=eig_args;
   args->cov_buffer=cov_buffer;
   args->eig_buffer=eig_buffer;
   args->fft_buffer=*fft_buffer;

   args->vec_len = vec_len;
   args->vec_num = vec_num;
   args->data = data;

   args->ifftFlag = ifftFlag;
   args->bits = bits;
   args->cov_mat_means=cov_mat_means;
   
   arm_rfft_fast_init_f32(&args->fft_buffer, vec_len);
   
   //change the vec leng only when its out of the fft, keep the length when its going into the matrix

   arm_mat_init_f32(args->input_buffer, vec_len, vec_num, data);
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
   arm_rfft_fast_f32(&args->fft_buffer,p + i*vec_len, pOut+ i*(vec_num/2+1), args->ifftFlag);
   }	
}

void cov(struct fft_pca_args* args)
//float32_t* A, float32_t* cov_mat, float32_t* means, uint32_t vec_len, uint32_t vec_num)
{
    float32_t* means = args->cov_mat_means;
    uint16_t vec_len = args->vec_len;
    uint16_t vec_num = args->vec_num;
    float32_t* cov_mat = args->cov_buffer;
    float32_t* inp_buf = (float32_t*) args->input_buffer;

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

   fft_obs_matrix(args);

   vec_len = (vec_len/2 + 1); // since data is real, vectors after fft are length n/2 + 1

   cov(args);
   eig_decomp(args->eig_args);
}

int32_t main(){
   while(1);
   struct fft_pca_args* args;
   fft_pca(args);
}

