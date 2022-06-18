#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "eig_vec_decomp.h"

void normalize(float* vec, uint32_t vec_len);
float l1_error(float* new_vec, float* old_vec, uint32_t vec_len);
void matrix_vec_mult(float* mat, uint32_t dim_size, float* vec, float* new_vec);


void print_matrix(float* mat, int columns, int rows) {
    for (int i = 0; i < columns*rows; i++) {
        printf("%f ", mat[i]);
        if ((i + 1) % columns == 0) {
            printf("\n");
        }
    }
    printf("\n");
}
void print_matrix_sci(float* mat, int columns, int rows) {
    for (int i = 0; i < columns*rows; i++) {
        printf("%e ", mat[i]);
        if ((i + 1) % columns == 0) {
            printf("\n");
        }
    }
    printf("\n");
}

void* alloc_eig_args(uint32_t dim_size, uint32_t eig_vec_num, uint32_t execs, float err_tol)
{
    struct eig_decomp_args* args = malloc(sizeof(struct eig_decomp_args));

    args->dim_size = dim_size;
    args->eig_vec_num = eig_vec_num;
    args->execs = execs;
    args->err_tol = err_tol;
    args->eig_vec = malloc(sizeof(float)*dim_size);
    args->s = malloc(sizeof(float)*dim_size);
    args->Sw = malloc(dim_size * sizeof(float));
    args->deflation_matrix = malloc(dim_size*dim_size*sizeof(float));

    for (int i = 0; i < dim_size; i++)
    {
        args->eig_vec[i] = 1;
    }

    return (void*)args;
}

void free_eig_args(struct eig_decomp_args* args) 
{
     free(args->s);
     free(args->eig_vec);
     free(args->Sw);
     free(args->deflation_matrix);
     free(args);
}

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

/*
 * inner_product
 *
 * Computes the inner (dot) product between vectors a and b, both of length
 * dim_size. Returns the scalar result.
 */
float inner_product(float *a, float *b, unsigned int dim_size) {
    float result = 0.;

    for(int k = 0; k < dim_size; k++) {
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
void outer_product(float *a, float *b, unsigned int dim_size, float *out) {
    for (int row = 0; row < dim_size; row++) {
        for (int col = 0; col < dim_size; col++) {
            out[row*dim_size+col] = a[row] * b[col];
        }
    }   
}

int power_iteration(float* matrix, struct eig_decomp_args* args)
{
    float* mat = matrix;
    float* eig_vec = args->eig_vec;
    float* s = args->s;
    uint32_t dim_size = args->dim_size;
    uint32_t execs = args->execs;
    float err_tol = args->err_tol;
    float err = 0;
    int i = 0;
    
    for (i = 0; i < dim_size; i++)
    {
        eig_vec[i] = 1;
    }

    for (i = 0; i < execs; i++)
    {
        matrix_vec_mult(mat, dim_size, eig_vec, s);
        normalize(s, dim_size);
        err = l1_error(s, eig_vec, dim_size);

        // swap the buffers for eigen vectors
        float* tmp = eig_vec;
        eig_vec = s;
        s = tmp;
        
        if (err < err_tol)
            break;
    }

    if (i == execs) {
        printf("COULDN'T CONVERGE, Error: %e\n", err);
    }

    args->eig_vec = eig_vec;
    args->s = s;

    return i;
}


void eig_decomp(float* matrix, float* eig_vecs, int* convergence, void* args) {
    int i,k;
    struct eig_decomp_args* e_args = (struct eig_decomp_args*) args;
    float* eig_vec = e_args->eig_vec;
    float* Sw = e_args->Sw;
    float* deflation_matrix = e_args->deflation_matrix;
    uint32_t nvecs = e_args->eig_vec_num;
    uint32_t dim_size = e_args->dim_size;
    
    for(k = 0; k < nvecs; k++) {

        // First, use power iteration to extract the dominant eigenvector of matrix.
        convergence[k] = power_iteration(matrix, e_args);
        
        // Copy eigenvector into output buffer for extracted eigenvectors
        memcpy(&eig_vecs[k*dim_size], eig_vec, dim_size * sizeof(float));

        // Warn if convergence is not reached
        if (convergence[k] == e_args->execs)
        {
            printf("EIGEN VECTOR: %d\n", k);
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
        matrix_vec_mult(matrix, dim_size, eig_vec, Sw);

        // inner product w'*S*w
        // verified correct result wTSw in Matlab
        float wTSw = inner_product(Sw, eig_vec, dim_size);

        // Element-wise multiplication.
        memcpy(Sw, eig_vec, dim_size * sizeof(float));
        for(i = 0; i < dim_size; i++) {
            Sw[i] *= wTSw;
        }

        // Finally, do the outer product between w and w'*S*w*w', which is stored in variable Sw
        outer_product(eig_vec, Sw, dim_size, deflation_matrix);

        // Last, subtract the deflation matrix from S.
        for(i = 0; i < dim_size*dim_size; i++){
            matrix[i] -= deflation_matrix[i];
        }
	
    }

}

/*
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

    float mat3[] = {5.3080155e+01,  4.1702272e+01,  4.5173333e+01,  3.7204400e+01,  2.7568818e+01,  1.2044329e+01,  3.0815801e+01,  4.1705880e+01, -2.6851760e+01,
                     4.1702272e+01,  5.2830289e+01,  3.8106160e+01,  4.4241230e+01,  3.6952886e+01,  2.3231513e+01,  4.0466756e+01,  4.3478852e+01, -1.3968003e+01,
                     4.5173333e+01,  3.8106160e+01,  5.9978970e+01,  3.0968549e+01,  1.6234313e+01, -8.2384427e-01,  2.4081191e+01,  3.5565437e+01, -3.9179657e+01,
                     3.7204400e+01,  4.4241230e+01,  3.0968549e+01,  5.4868431e+01,  4.6475526e+01,  3.5481616e+01,  4.6860827e+01,  4.6418166e+01, -2.2520432e-01,
                     2.7568818e+01,  3.6952886e+01,  1.6234313e+01,  4.6475526e+01,  6.1248984e+01,  4.5543887e+01,  4.9276593e+01,  4.5075905e+01,  1.8513892e+01,
                     1.2044329e+01,  2.3231513e+01, -8.2384427e-01,  3.5481616e+01,  4.5543887e+01,  5.1558562e+01,  4.0096637e+01,  2.9939223e+01,  2.8770362e+01,
                     3.0815801e+01,  4.0466756e+01,  2.4081191e+01,  4.6860827e+01,  4.9276593e+01,  4.0096637e+01,  5.7204552e+01,  4.6293188e+01,  1.0893214e+01,
                     4.1705880e+01,  4.3478852e+01,  3.5565437e+01,  4.6418166e+01,  4.5075905e+01,  2.9939223e+01,  4.6293188e+01,  5.4687541e+01, -6.8967558e+00,
                    -2.6851760e+01, -1.3968003e+01, -3.9179657e+01, -2.2520432e-01,  1.8513892e+01,  2.8770362e+01,  1.0893214e+01, -6.8967558e+00,  5.9540503e+01};

    uint32_t dim_size = 9;
    uint32_t execs = 100000;
    float err_tol = 0.001;

    struct eig_decomp_args* args = alloc_eig_args(dim_size,execs, err_tol);

    power_iteration(mat1, args);

    float *vecs = malloc(dim_size*dim_size*sizeof(float));
    eig_decomp(mat3, vecs, 3, args); */
    /* should be close to V1 = {-0.1380, 0.2371, 0.6335,  0.1986,  0.0152, 0.3878, 0.1201,  0.5648} */
    /* should be close to V2 = {-0.1625, 0.0062, 0.3589, -0.4096, -0.0023, 0.5057, 0.0039, -0.6489} */
    /* should be close to V2 =  */
    /* should be close to V2 =  */
//    print_matrix(args->eig_vec, dim_size, 1);
//    print_matrix(vecs, dim_size, 3);

//    args->targ_mat = mat2;

//    eig_decomp(args);

    /* should be close to {-0.2137, -0.3353, -0.0551, -0.6872, 0.3163, -0.3101, 0.1152, 0.3964} */
//    print_matrix(args->eig_vec, dim_size, 1);

//    return 1;
//}

