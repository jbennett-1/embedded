#include "arm_math.h"
#include <math.h>
#define ARMCM7

struct eig_decomp_args {
    float32_t*   s;         // mem needed to calc error
    arm_matrix_instance_f32 * targ_mat;
    uint32_t dim_size;
    float32_t*   eig_vec;  
    uint32_t execs;    // maxmimum loop executions before giving up on convergence
    float32_t    err_tol;  // error tolerated in the eigenvector
};

void normalize(float32_t* vec, uint32_t vec_len)
{   
    float32_t norm = 0;
    
    arm_dot_prod_f32(vec, vec, vec_len, &norm);
    
    norm = sqrtf(norm);
    
    arm_scale_f32(vec, 1/norm, vec, vec_len);
}

float32_t l1_error(float32_t* new_vec, float32_t * old_vec, uint32_t vec_len)
{   
    float32_t error = 0;
    
    for (int i = 0; i < vec_len; i++)
    {   
        error += fabsf(new_vec[i] - old_vec[i]);
    }
    
    return error;
}

void eig_decomp(struct eig_decomp_args* args)
{   
    arm_matrix_instance_f32* mat = args->targ_mat;
    float32_t* eig_vec = args->eig_vec;
    float32_t* s = args->s;
    uint32_t dim_size = args->dim_size;
    uint32_t execs = args->execs;
    float32_t err_tol = args->err_tol;

    for (int i = 0; i < execs; i++)
    {
        arm_mat_vec_mult_f32(mat, eig_vec, s);
        normalize(s, dim_size);
        float32_t err = l1_error(s, eig_vec, dim_size);

        // swap the buffers for eigen vectors
        float32_t* tmp = eig_vec;
        eig_vec = s;
        s = tmp;

        if (err < err_tol)
            break;
    }

    args->eig_vec = eig_vec;
}

int32_t main(){
	while(1);
}
