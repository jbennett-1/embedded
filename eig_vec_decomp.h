#ifndef EIG_VEC_DECOMP_H
#define EIG_VEC_DECOMP_H

struct eig_decomp_args {
    float*   s;         // mem needed to calc error
    float*   Sw;
    float*   deflation_matrix;
    float*   eig_vec;
    uint32_t dim_size;
    uint32_t eig_vec_num;
    uint32_t execs;    // maxmimum loop executions before giving up on convergence
    float    err_tol;  // error tolerated in the eigenvector
};

void free_eig_args(struct eig_decomp_args* args);
void* alloc_eig_args(uint32_t dim_size, uint32_t eig_vec_num, uint32_t execs, float err_tol);
void eig_decomp(float* matrix, float* eig_vecs, int* convergence, void* args);

#endif
