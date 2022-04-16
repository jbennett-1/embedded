#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

/*
    Python code in case you wanna look at that

    import numpy as np
    from numpy.linalg import norm

    def eig_decomp(X, vec_len, vec_num, c, tol):
    
    r = np.ones(vec_len)
    r = r / norm(r, 2)

    for i in range(c):
        s = X.dot(r)
        s = s / norm(s, 2)
        error = norm(s - r, 1)
        r = s
        if error < tol:
            break
    return r
*/

struct eig_decomp_args {
    float*   s;         // mem needed to calc error
    float*   targ_mat;
    uint32_t dim_size;
    float*   eig_vec;
    uint32_t execs;    // maxmimum loop executions before giving up on convergence
    float    err_tol;  // error tolerated in the eigenvector
};

void normalize(float* vec, uint32_t vec_len)
{
    float norm = 0;

    for (int i = 0; i < vec_len; i++)
    {
        norm += vec[i]*vec[i];
    }

    norm = sqrtf(norm);

    for (int i = 0; i < vec_len; i++)
    {
        vec[i] = vec[i] / norm;
    }
}

float l1_error(float* new_vec, float* old_vec, uint32_t vec_len)
{
    float error = 0;

    for (int i = 0; i < vec_len; i++)
    {
        error += fabsf(new_vec[i] - old_vec[i]);
    }

    return error;
}

void matrix_vec_mult(float* mat, uint32_t dim_size, float* vec, float* new_vec)
{
    for (int i = 0; i < dim_size; i++)
    {
        new_vec[i] = 0;
        for (int j = 0; j < dim_size; j++)
        {
            new_vec[i] += mat[i*dim_size + j] * vec[j];
        }
    }
}

void eig_decomp(struct eig_decomp_args* args)
{
    float* mat = args->targ_mat;
    float* eig_vec = args->eig_vec;
    float* s = args->s;
    uint32_t dim_size = args->dim_size;
    uint32_t execs = args->execs;
    float err_tol = args->err_tol;

    for (int i = 0; i < execs; i++)
    {
        matrix_vec_mult(mat, dim_size, eig_vec, s);
        normalize(s, dim_size);
        float err = l1_error(s, eig_vec, dim_size);

        // swap the buffers for eigen vectors
        float* tmp = eig_vec;
        eig_vec = s;
        s = tmp;
        
        if (err < err_tol)
            break;
    }

    args->eig_vec = eig_vec;
}

struct eig_decomp_args* alloc_args(float* mat, uint32_t dim_size, uint32_t execs, float err_tol)
{
    struct eig_decomp_args* args = malloc(sizeof(struct eig_decomp_args));

    args->targ_mat = mat;
    args->dim_size = dim_size;
    args->execs = execs;
    args->err_tol = err_tol;
    args->eig_vec = malloc(sizeof(float)*dim_size);
    args->s = malloc(sizeof(float)*dim_size);

    for (int i = 0; i < dim_size; i++)
    {
        args->eig_vec[i] = 1;
    }

    normalize(args->eig_vec, dim_size);

    return args;
}

void print_matrix(float* mat, int columns, int rows) {
    for (int i = 0; i < columns*rows; i++) {
        printf("%f ", mat[i]);
        if ((i + 1) % columns == 0) {
            printf("\n");
        }
    }
    printf("\n");
}

// Test program
int main() {

     float mat1[] = {1.1763691e-02, -2.0218568e-02, -5.4010831e-02, -1.6934412e-02, -1.2958426e-03, -3.3064850e-02, -1.0242928e-02, -4.8154078e-02,
                    -2.0218568e-02,  3.4750186e-02,  9.2829846e-02,  2.9105622e-02, 2.2271990e-03,  5.6829434e-02,  1.7604792e-02,  8.2763694e-02,
                    -5.4010831e-02,  9.2829846e-02,  2.4798085e-01,  7.7751249e-02, 5.9496239e-03,  1.5181120e-01,  4.7028527e-02,  2.2109064e-01,
                    -1.6934412e-02,  2.9105622e-02,  7.7751249e-02,  2.4377920e-02, 1.8654292e-03,  4.7598477e-02,  1.4745199e-02,  6.9320172e-02,
                    -1.2958426e-03,  2.2271990e-03,  5.9496239e-03,  1.8654292e-03, 1.4274500e-04,  3.6422957e-03,  1.1283213e-03,  5.3044669e-03,
                    -3.3064850e-02,  5.6829434e-02,  1.5181120e-01,  4.7598477e-02, 3.6422957e-03,  9.2937179e-02,  2.8790358e-02,  1.3534930e-01,
                    -1.0242928e-02,  1.7604792e-02,  4.7028527e-02,  1.4745199e-02, 1.1283213e-03,  2.8790358e-02,  8.9187631e-03,  4.1928913e-02,
                    -4.8154078e-02,  8.2763694e-02,  2.2109064e-01,  6.9320172e-02, 5.3044669e-03,  1.3534930e-01,  4.1928913e-02,  1.9711632e-01};

    float mat2[] = {0.01769675, 0.02776701, 0.00456199, 0.05690041, -0.02619395, 0.02567673, -0.00953805, -0.03282689,
                    0.02776701, 0.04356773, 0.00715798, 0.08927937, -0.04109953, 0.04028798, -0.01496565, -0.0515069,
                    0.00456199, 0.00715798, 0.00117602, 0.01466819, -0.00675246, 0.00661913, -0.00245879, -0.00846235,
                    0.05690041, 0.08927937, 0.01466819, 0.18295205, -0.0842215,  0.08255848, -0.03066773, -0.10554838,
                    -0.02619395, -0.04109953, -0.00675246, -0.0842215, 0.03877115, -0.03800558, 0.01411781,  0.04858892,
                    0.02567673, 0.04028798, 0.00661913, 0.08255848, -0.03800558, 0.03725513, -0.01383904, -0.04762949,
                    -0.00953805, -0.01496565, -0.00245879, -0.03066773, 0.01411781, -0.01383904, 0.00514074,  0.01769277,
                    -0.03282689, -0.0515069,  -0.00846235, -0.10554838,  0.04858892, -0.04762949, 0.01769277,  0.06089279};

    uint32_t dim_size = 8;
    uint32_t execs = 1000;
    float err_tol = 0.001;

    struct eig_decomp_args* args = alloc_args(mat1, dim_size, execs, err_tol);

    eig_decomp(args);

    /* should be close to {-0.1380, 0.2371, 0.6335, 0.1986, 0.0152, 0.3878, 0.1201, 0.5648} */
    print_matrix(args->eig_vec, dim_size, 1);

    args->targ_mat = mat2;

    eig_decomp(args);
i
    /* should be close to {-0.2137, -0.3353, -0.0551, -0.6872, 0.3163, -0.3101, 0.1152, 0.3964} */
    print_matrix(args->eig_vec, dim_size, 1);

    return 1;
}
