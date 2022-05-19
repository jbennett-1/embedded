#include "arm_math.h"
#include "eig_vec_decomp.h"
#define ARMCM7_DP

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
//	y  = y * ( threehalfs - ( x2 * y * y ) );   // 2nd iteration, this can be removed

	return y;
}

void normalize(float32_t* vec, uint32_t vec_len)
{   
    float32_t norm = 0;
    
    for(int i = 0; i < vec_len; i++)
    {
	arm_dot_prod_f32(vec, vec, vec_len, &norm);
    }

    norm = Q_rsqrt(norm);

    for(int i = 0; i < vec_len; i++)
    {    
	arm_scale_f32(vec, 1/norm, vec, vec_len);
    }

}

float32_t l1_error(float32_t* new_vec, float32_t * old_vec, uint32_t vec_len)
{   
    float32_t error = 0;
    
    for (int i = 0; i < vec_len; i++)
    {   
	float32_t diff = new_vec[i] - old_vec[i];
	error += (diff < 0) ? (diff*-1) : (diff); 
	}
    
    return error;
}

void eig_decomp(arm_matrix_instance_f32* matrix, struct eig_decomp_args* args)
{   
    args->targ_mat=matrix;
    float32_t* eig_vec = args->eig_vec;
    float32_t* s = args->s;
    uint32_t dim_size = args->dim_size;
    uint32_t execs = args->execs;
    float32_t err_tol = args->err_tol;
    uint64_t i;
    float32_t err;

    for(int16_t j = 0; j < dim_size; j++){
	eig_vec[j]=1;
    }

    for (i = 0; i < execs; i++)
    {
        arm_mat_vec_mult_f32(matrix, eig_vec, s);
        normalize(s, dim_size);
        err = l1_error(s, eig_vec, dim_size);

        // swap the buffers for eigen vectors
        float32_t* tmp = eig_vec;
        eig_vec = s;
        s = tmp;

        if (err < err_tol){
            break;
	}
    }
    args->eig_vec = eig_vec;
    args->s=s;
}

