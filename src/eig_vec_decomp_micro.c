#include "arm_math.h"
#include "system.h"
#include "eig_vec_decomp.h"
#define ARMCM4_FP
#include <math.h>

void normalize(float32_t* vec, uint32_t vec_len);
float32_t l1_error(float32_t* new_vec, float32_t* old_vec, uint32_t vec_len);
void matrix_vec_mult(float32_t* mat, uint32_t dim_size, float32_t* vec, float32_t* new_vec);

void bit_gen(float32_t* eig_vec,struct eig_decomp_args* e_args){
    uint32_t mean, max, min, sum, bit_adder=0;
    uint16_t i, k, t;
    uint32_t dim_size = e_args->dim_size;    
    
    for(i = 0;i < dim_size * e_args->eig_vec_num;i++) {
	sum += eig_vec[i];
        mean = sum / (dim_size*4);
    }

    for(k = 0; k < dim_size*4; k++) {
	if(eig_vec[i] > max) { //max
	    max = eig_vec[k];
    	}

	if(eig_vec[i] < min) { //min
	    min = eig_vec[k];
	}
     }

    for(t = 0;t < dim_size*4; t++) {
	bit_adder = eig_vec[t];
        bit_adder |= (mean>0) ? 1 : 0;
        bit_adder <<=1;
        bit_adder <<=0;

        if(bit_adder > (max + mean) / 2) {
	    bit_adder |=  1 << 1;
    	} else {
	    bit_adder |= 0 << 1;
	}

	if(bit_adder |= 0 << 0) {
	    if(bit_adder > (max + mean) /2) {
	 	bit_adder |= 1 << 1;
	    } else {
		bit_adder |= 1 << 1;
	    }
	}
    }   
}

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

void matrix_vec_mult(float32_t* mat, uint32_t dim_size, float32_t* vec, float32_t* new_vec)
{
    for (int32_t i = 0; i < dim_size; i++)
    {
        new_vec[i] = 0;
        for (int32_t j = 0; j < dim_size; j++)
        {
            new_vec[i] += mat[i*dim_size + j] * vec[j];
        }
    }
}

float32_t l1_error(float32_t* new_vec, float32_t* old_vec, uint32_t vec_len)
{
    float32_t error = 0;

    for (int32_t i = 0; i < vec_len; i++)
    {
    error = (new_vec[i]-old_vec[i] < 0) ? -(new_vec[i]-old_vec[i]) : (new_vec[i] - old_vec[i]);
//    error += fabsf(new_vec[i] - old_vec[i]);
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

int power_iteration(arm_matrix_instance_f32* matrix, struct eig_decomp_args* e_args)
{
    float32_t* eig_vec = e_args->eig_vec;
    float32_t* s = e_args->s;
    uint32_t dim_size = e_args->dim_size;
    uint32_t execs = e_args->execs;
    float32_t err_tol = e_args->err_tol;
    float32_t err = 0;
    uint32_t i = 0;
   
    for (int i = 0; i < dim_size; i++)
    {
        eig_vec[i] = 1;
    }
 
    
    for (i = 0; i < execs; i++)
    {
        matrix_vec_mult(matrix->pData, dim_size, eig_vec, e_args->s);
	normalize(s, dim_size);
        err = l1_error(s, eig_vec, dim_size);

        // swap the buffers for eigen vectors
        float32_t* tmp = eig_vec;
        eig_vec = s;
        s = tmp;
        
        if (err < err_tol)
            break;
    }

    eig_vec = eig_vec;
    s = s;

    return i;
}


void eig_decomp(arm_matrix_instance_f32* matrix, struct eig_decomp_args* e_args) {
    float32_t* eig_vec = e_args->eig_vec;
    float32_t* Sw = e_args->Sw;
    float32_t* deflation_matrix = e_args->deflation_matrix;
    int32_t i,k;
    uint32_t nvecs = e_args->eig_vec_num;
    uint32_t dim_size = e_args->dim_size;
    
    for(k = 0; k < nvecs; k++) {

        // First, use power iteration to extract the dominant eigenvector of matrix.
        power_iteration(matrix, e_args);
        // Copy eigenvector into output buffer for extracted eigenvectors
        for(int s =0; s<dim_size;s++) {
		eig_vec[k*dim_size+s]=eig_vec[s];
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
        matrix_vec_mult(matrix->pData, dim_size, eig_vec, Sw);

        // inner product w'*S*w
        // verified correct result wTSw in Matlab
        float32_t wTSw = inner_product(Sw, eig_vec, dim_size);

        // Element-wise multiplication.
        for(int32_t t=0;t<dim_size;t++){
		Sw[t]=eig_vec[t];
	}
        
	for(i = 0; i < dim_size; i++) {
            Sw[i] *= wTSw;
        }

        // Finally, do the outer product between w and w'*S*w*w', which is stored in variable Sw
        outer_product(eig_vec, Sw, dim_size, deflation_matrix);

        // Last, subtract the deflation matrix from S.
        for(i = 0; i < dim_size*dim_size; i++){
            matrix->pData[i]-= deflation_matrix[i];
        }
	
    }

}
