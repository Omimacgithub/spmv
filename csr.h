#ifndef CSR_MATRIX_H
#define CSR_MATRIX_H
#include <stddef.h>

struct csr{
    size_t *row_offsets;
    double *values;
    size_t *column_indices;
};

typedef struct csr csr;


csr intoCSR(unsigned int size, double *mat);

#endif //CSR_MATRIX_H
