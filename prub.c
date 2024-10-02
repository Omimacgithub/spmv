#include "csr.h"
#include <stdio.h>

void main(){
  double mat[] = {1.0, 0.0, 2.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 5.0, 0.0, 0.0, 0.0, 0.0, 6.0, 7.0, 8.0, 0.0};
  csr m = intoCSR(5, mat);
  int i;
  printf("m.row_offsets:\n {");
  for (i=0; i < 5; i++){
	printf("%d, ", m.row_offsets[i]);
  }
  printf("}\n");
  int valsize = m.row_offsets[i];
  printf("m.column_indices:\n {");
  for (i=0; i < valsize; i++){
	printf("%d, ", m.column_indices[i]);
  }
  printf("}\n");
  printf("m.values:\n {");
  for (i=0; i < valsize; i++){
	printf("%f, ", m.values[i]);
  }
  printf("}\n");
 

}
