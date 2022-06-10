#ifndef EIG_VEC_DECOMP_H
#define EIG_VEC_DECOMP_H
#include "arm_math.h"

struct eig_decomp_args {
    float* s;         // mem needed to calc error
    uint32_t dim_size;
    float32_t* eig_vec;
    uint32_t execs;    // maxmimum loop executions before giving up on convergence
    float32_t err_tol;  // error tolerated in the eigenvector
};

void eig_decomp(arm_matrix_instance_f32* targ_mat, struct eig_decomp_args* eig_args);

#endif
