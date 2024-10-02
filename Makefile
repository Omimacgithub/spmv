spmv: link

	./spmv

link: compile

	gcc my_dense.o my_sparse.o timer.o spmv.o -lopenblas -lm -o spmv

compile: 

	make my_dense.o my_sparse.o timer.o spmv.o
