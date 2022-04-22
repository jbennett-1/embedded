#include <stdint.h>
#include "arm_math.h"
	
struct fft_pca_args{
	void fft_args;
	arm_matrix_instance_f32* input_buffer;
	float32_t* eig_buffer;
	arm_matrix_instance_f32* cov_buffer;
	float32_t* output_buffer;
}

//dynamic mem management//calloc() but write your own version//alloc at compile time//constants in ld for mem 
//need 2 buffers, data, cov, and 2 others, eigen vector buffer + store the product in output buffer
/* matrix is stored as an array, contigious memory in stack, need pointer array (for the buffers of elements x1-x2)
//using headerfiles definitions of arm_matrix instance etc etc
for loop to go til end of vector...
use rfft from the transformFunctions will produce spaces between the buffer,, get rid of those,, it goes to cov last */

void initialize(struct fft_pca_args* args)
{
	arm_matrix_instance_f32* mat_buf = args->input_buffer;
	uint8_t ifftFlag;
	arm_rfft_fast_instance_f32 *fft_buf = args->input_buffer;

	uint16_t vec_len = args->vec_len;

	arm_rfft_fast_init_f32(*fft_buf, 
	arm_rfft_fast_f32(*fft_buf, args->input_buffer, args->output_buffer, ifftFlag);
}
void fft_pca(uint8_t* bits, struct fft_pca_args* args)
{
	// ???
        // TODO: put data_size in args
//	arm_mat_init_f32(arm_matrix_instance_f32 *mat_buf, uint16_t vec_len, uint16_t vec_num, float32_t* data_size);
	arm_rfft_fast_instance_f32 *fft_buf;
	uint8_t ifftFlag;

	arm_rfft_fast_f32(*fft_buf, args->input_buffer, args->output_buffer, ifftFlag);

	//arm_rfft_fast_instance_f32* fft_buf = inputs->src_rfft_mat;
	//arm_matrix_instance_f32* args->src_mat=input_buffer;
	//float32_t* input_buffer = args->src_mat;	
	
	float32_t* cov_buf = args->cov_mat;
	float32_t* cov_mat_means = args->cov_mat_means;
	float32_t* data_size;	

	uint16_t vec_len = args->vec_len;
	vec_len = (vec_len/2 + 1); // since data is real, vectors after fft are length n/2 + 1
	uint16_t vec_num = args->vec_num;


        arm_mat_init_f32(*mat_buf, vec_len, vec_num, data_size);
	//input_buffer is type arm_matrix_instance // vec_len is uint32 // vec_num is uint32
	fft_obs_matrix(*mat_buf, vec_len, vec_num);	
	//input_buffer is arm_matrix_instance // cov_buf is float32_t // cov_mat_means is float 32 // vec num is uint32
	cov(args->input_buffer, args->cov_buffer, cov_mat_means, *mat_buf, vec_num);
}

/*
  void arm_rfft_fast_f32(
        const arm_rfft_fast_instance_f32 * S,
        float32_t * p, float32_t * pOut,
        uint8_t ifftFlag);

void arm_mat_init_f32(
        arm_matrix_instance_f32 * S,
        uint16_t nRows,
        uint16_t nColumns,
        float32_t * pData);


*/
/*	ssyevx_("Vectors", "Indices", "Upper", (integer*) &n, (real*) cov_mat, (integer*) &lda, (real*) 0, (real*) 0,
		(integer*) &upper_limit, (integer*) &lower_limit, (real*) &absol, (integer*) &eig_val_found, (real*) &eig_val,
		(real*) princ_comp, (integer*) &ldz, (real*) work, (integer*) &lwork, (integer*) iwork, (integer*) &ifail,
		(integer*) &info);
*/



void fft_obs_matrix(float32_t* mat_buf, struct fft_pca_args* args)
{
	arm_matrix_instance_f32* args->src_mat = mat_buf;
	
	for (int i = 0; i < vec_num; i++)
	{
	/*  Since the data is real, the output should be length n/2 + 1. So,
            placing the output as (buffer + i*(vec_len/2 + 1)) allows us to 
            condense the matrix without leaving any extra space because of
            how real ffts work. */
        arm_rfft_fast_f32((mat_buf + i*vec_len), (mat_buf + i*(vec_len/2 + 1)), args);
    }
}

void cov(float32_t* A, float32_t* cov_mat, float32_t* means, uint32_t vec_len, uint32_t vec_num)
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
