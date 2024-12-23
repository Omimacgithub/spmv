TARGETS = gsl mkl
CC = gcc
ICC = icx
OFLAGS = -march=native -g3 -O3
WARNFLAGS = -Wall -Wextra
CFLAGS =  $(OFLAGS) $(WARNFLAGS) -D_GSL_
IFLAGS =  $(OFLAGS) $(WARNFLAGS) -D_MKL_ #-I"${MKLROOT}/include" -diag-disable=10441
ADVFLAGS = -Rpass-analysis=loop-vectorize -Rpass=loop-vectorize -Rpass-missed=loop-vectorize -qopt-report=3 -vec -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -liomp5 -lpthread 
#ADVFLAGS = -Rpass-analysis=loop-vectorize -Rpass=loop-vectorize -Rpass-missed=loop-vectorize -qopt-report=3 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -liomp5
LDLIBS = -lm -lgsl -lgslcblas
ILDLIBS = -lm -lmkl

SRCS = timer.c my_dense.c my_csr.c my_coo.c my_csc.c spmv.c
IOBJS = $(SRCS:.c=.o)
OBJS = $(SRCS:.c=_g.o)

all: $(TARGETS)

gsl: $(OBJS) 	
	$(CC) $(CFLAGS) $^ $(LDLIBS) -o $@
mkl: $(IOBJS)
	$(ICC) $(IFLAGS) $(ADVFLAGS) $^ $(ILDLIBS) -o $@

$(IOBJS): $(SRCS)
	$(ICC) $(IFLAGS) $(ADVFLAGS) $^ -c

%_g.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

cleanall:
	$(RM) $(OBJS) $(IOBJS) $(TARGETS) *~

