#include "csr.h"
#include <stdlib.h>

struct Flag {
    unsigned int flag : 1;  // 1-bit field
};


csr intoCSR(unsigned int size, double *mat){
	unsigned int i, j=0, z=0, msize = size*size;
	csr m;
	struct Flag f;

	f.flag=1;

	//m.ros = (size_t) size+1;
	m.values = (double *) malloc(size * size * sizeof(double));
	m.row_offsets = (size_t *) malloc( (size+1) * sizeof(size_t));
	m.column_indices = (size_t *) malloc( size * size * sizeof(size_t));

	for(i=0; i < msize; i++){
		if (i%size == 0){
			//We're in a new row of the matrix
			if(f.flag == 1 && i != 0){
				//There were no values on the entire row
				m.row_offsets[z] = m.row_offsets[z-1];
				++z;
			}
			f.flag = 1;
		}
		if(mat[i] != 0){
			m.values[j] = mat[i];
			m.column_indices[j] = i%size;
			if (f.flag == 1){
				//This value is the first ocurrence of the row
				f.flag = 0;
				m.row_offsets[z] = j;
				++z;
			}
			++j;
		}
	}


	//Resize vectors
	m.values = (double *) realloc(m.values, j * sizeof(double));
	m.column_indices = (size_t *) realloc(m.column_indices, j * sizeof(size_t));
	m.row_offsets[z] = j;

	return m;
}
