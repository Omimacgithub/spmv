#include "csr.h"
#include <stdlib.h>

struct Flag {
    unsigned int flag : 1;  // 1-bit field
};


csr intoCSR(unsigned int size, double *mat){
	int i, j=0, z=0, msize = size*size;
	csr m;
	struct Flag f;

	f.flag=1;

	m.values = (double *) malloc(size * size * sizeof(double));
	m.row_offsets = (int *) malloc( (size+1) * sizeof(int));
	m.column_indices = (int *) malloc( size * size * sizeof(int));

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
		if(mat[i] != 0.0){
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
	m.column_indices = (int *) realloc(m.column_indices, j * sizeof(int));
	//The last elemente is the size of the values vector
	m.row_offsets[z] = j;

	return m;
}
