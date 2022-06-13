#include "arm_math.h"
#include "system.h"
#include "eig_vec_decomp.h"
#define ARMCM4_FP
#include <math.h>
/*
float Q_rsqrt( float number )
{
        long i;
        float x2, y;
        const float threehalfs = 1.5F;

        x2 = number * 0.5F;
        y  = number;
        i  = * ( long * ) &y;                       // evil floating point bit level hacking
        i  = 0x5f3759df - ( i >> 1 );               // what the fuck? 
        y  = * ( float * ) &i;
        y  = y * ( threehalfs - ( x2 * y * y ) );   // 1st iteration
//      y  = y * ( threehalfs - ( x2 * y * y ) );   // 2nd iteration, this can be removed
        return y;

}
*/

void normalize(float32_t* vec, uint32_t vec_len)
{
    float32_t norm = 0;
    
    for(int32_t i = 0; i < vec_len; i++)
    {
	arm_dot_prod_f32(vec, vec, vec_len, &norm);
    }

    norm = sqrt(norm);

    for(int32_t i = 0; i < vec_len; i++)
    {    
	arm_scale_f32(vec, 1/norm, vec, vec_len);
    }

}

float32_t l1_error(float32_t* new_vec, float32_t * old_vec, uint32_t vec_len)
{   
    float32_t error = 0;
    
    for (int32_t i = 0; i < vec_len; i++)
    {   
	float32_t diff = new_vec[i] - old_vec[i];
	error += (diff < 0) ? (diff*-1) : (diff); 
	}
    
    return error;
}

void eig_decomp(arm_matrix_instance_f32* targ_mat, struct eig_decomp_args* eig_args)
{   
    uint64_t i;
    float32_t err;

    for(int16_t j = 0; j < eig_args->dim_size; j++){
	eig_args->eig_vec[j]=1;
    }

    for (i = 0; i < eig_args->execs; i++)
    {
        arm_mat_vec_mult_f32(targ_mat, eig_args->eig_vec, eig_args->s);
        normalize(eig_args->s, eig_args->dim_size);
        err = l1_error(eig_args->s, eig_args->eig_vec, eig_args->dim_size);

        // swap the buffers for eigen vectors
        float32_t* tmp = eig_args->eig_vec;
        eig_args->eig_vec = eig_args->s;
        eig_args->s = tmp;

        if (err < eig_args->err_tol){
            break;
	}
    }
    eig_args->eig_vec = eig_args->eig_vec;
    eig_args->s=eig_args->s;
}

