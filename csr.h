#ifndef CSR_MATRIX_H
#define CSR_MATRIX_H

struct csr{
    int *row_offsets;
    double *values;
    int *column_indices;
};

typedef struct csr csr;


csr intoCSR(unsigned int size, double *mat);

#endif //CSR_MATRIX_H
