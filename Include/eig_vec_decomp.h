#ifndef EIG_VEC_DECOMP_H
#define EIG_VEC_DECOMP_H
#include "arm_math.h"
#include "system.h"
#include "eig_vec_decomp.h"

struct eig_decomp_args {
    float32_t*   s;         // mem needed to calc error
    float32_t*   Sw;
    float32_t*   deflation_matrix;
    float32_t*   eig_vec;
    uint32_t dim_size;
    uint32_t eig_vec_num;
    uint32_t execs;    // maxmimum loop executions before giving up on convergence
    float32_t    err_tol;  // error tolerated in the eigenvector
};

void eig_decomp(arm_matrix_instance_f32* targ_mat, struct eig_decomp_args* e_args);

#endif
