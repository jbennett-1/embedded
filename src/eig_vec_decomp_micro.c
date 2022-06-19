#include "arm_math.h"
#include "system.h"
#include "eig_vec_decomp.h"
#define ARMCM4_FP
#include <math.h>

void normalize(float32_t* vec, uint32_t vec_len);
float l1_error(float32_t* new_vec, float32_t* old_vec, uint32_t vec_len);
void matrix_vec_mult(float32_t* mat, uint32_t dim_size, float32_t* vec, float32_t* new_vec);

void normalize(float32_t* vec, uint32_t vec_len)
{
    float32_t norm = 0;

    for (int32_t i = 0; i < vec_len; i++)
    {
        norm += vec[i]*vec[i];
    }

    norm = sqrtf(norm);

    for (int32_t i = 0; i < vec_len; i++)
    {
        vec[i] = vec[i] / norm;
    }
}

float32_t l1_error(float32_t* new_vec, float32_t* old_vec, uint32_t vec_len)
{
    float32_t error = 0;

    for (int32_t i = 0; i < vec_len; i++)
    {
        error += fabsf(new_vec[i] - old_vec[i]);
    }

    return error;
}

/*
 * inner_product
 *
 * Computes the inner (dot) product between vectors a and b, both of length
 * dim_size. Returns the scalar result.
 */
float32_t inner_product(float32_t *a, float32_t *b, uint32_t dim_size) {
    float32_t result = 0.;

    for(int32_t k = 0; k < dim_size; k++) {
        result += a[k] * b[k];
    }
    return result;
}

/*
 * outer_product
 *
 * Computes the outer product between matricies a and b:
 *
 *    a  *  b  =  out
 *  (nx1) (1xn)  (nxn)
 *
 */
void outer_product(float32_t *a, float32_t *b, uint32_t dim_size, float32_t *out) {
    for (uint32_t row = 0; row < dim_size; row++) {
        for (uint32_t col = 0; col < dim_size; col++) {
            out[row*dim_size+col] = a[row] * b[col];
        }
    }   
}

int power_iteration(arm_matrix_instance_f32* mat, struct eig_decomp_args* args)
{
    float32_t* eig_vec = args->eig_vec;
    float32_t* s = args->s;
    uint32_t dim_size = args->dim_size;
    uint32_t execs = args->execs;
    float32_t err_tol = args->err_tol;
    float32_t err = 0;
    uint32_t i = 0;
    
    for(int32_t j = 0; j < dim_size; j++){
        eig_vec[j]=1;
    }

    
    for (i = 0; i < execs; i++)
    {
        arm_mat_vec_mult_f32(mat, eig_vec, s);
        normalize(s, dim_size);
        err = l1_error(s, eig_vec, dim_size);

        // swap the buffers for eigen vectors
        float32_t* tmp = eig_vec;
        eig_vec = s;
        s = tmp;
        
        if (err < err_tol)
            break;
    }

    args->eig_vec = eig_vec;
    args->s = s;

    return i;
}


void eig_decomp(arm_matrix_instance_f32* matrix, struct eig_decomp_args* e_args) {
    int32_t i,k;
    uint32_t nvecs = e_args->eig_vec_num;

    
    for(k = 0; k < nvecs; k++) {

        // First, use power iteration to extract the dominant eigenvector of matrix.
        power_iteration(matrix, e_args);
        // Copy eigenvector into output buffer for extracted eigenvectors
        for(int s =0; s<e_args->dim_size;s++) {
		e_args->eig_vec[k*e_args->dim_size+s]=e_args->eig_vec[s];
	}

        // Deflate the dominant eigenvector from the matrid
        // S = S - w*w'*S*w*w'
        // S = S - w*(w'*(S*w)*w'
        //                mat-vec-mult
        //           |  inner prod    |
        //           |  elem-wise mult     |
        //         |   outer product      |
        // matrix vector multiplication S*w
        // verified correct result (Sw) in Matlab
        arm_mat_vec_mult_f32(matrix, e_args->eig_vec, e_args->Sw);

        // inner product w'*S*w
        // verified correct result wTSw in Matlab
        float32_t wTSw = inner_product(e_args->Sw, e_args->eig_vec, e_args->dim_size);

        // Element-wise multiplication.
        for(int32_t t=0;t<e_args->dim_size;t++){
		e_args->Sw[t]=e_args->eig_vec[t];
	}
	// memcpy(Sw, eig_vec, dim_size * sizeof(float));
        
	for(i = 0; i < e_args->dim_size; i++) {
            e_args->Sw[i] *= wTSw;
        }

        // Finally, do the outer product between w and w'*S*w*w', which is stored in variable Sw
        outer_product(e_args->eig_vec, e_args->Sw, e_args->dim_size, e_args->deflation_matrix);

        // Last, subtract the deflation matrix from S.
        for(i = 0; i < e_args->dim_size*e_args->dim_size; i++){
            matrix->pData[i]-= e_args->deflation_matrix[i];
        }
	
    }

}
