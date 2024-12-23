TARGETS = gsl mkl
CC = gcc
ICC = icc
OFLAGS = -O0
WARNFLAGS = -Wall -Wextra
CFLAGS = -g3 $(OFLAGS) $(WARNFLAGS) -D_GSL_
IFLAGS = -g3 $(OFLAGS) $(WARNFLAGS) -D_MKL_ -I"${MKLROOT}/include" -diag-disable=10441
ADVFLAGS = -qopt-report=5 -vec -simd -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -liomp5 -lpthread
LDLIBS = -lm -lgsl -lgslcblas
ILDLIBS = -lm -lmkl
#xCORE-AVX2 xCOMMON-AVX512 <- advisor recomienda esto 

SRCS = timer.c my_dense.c my_csr.c my_coo.c my_csc.c spmv.c
IOBJS = $(SRCS:.c=.o)
OBJS = $(SRCS:.c=_g.o)

all: $(TARGETS)

gsl: $(OBJS) 	
	$(CC) $(CFLAGS) $^ $(LDLIBS) -o $@
mkl: $(IOBJS)
	$(ICC) $(IFLAGS) $(ADVFLAGS) $^ $(ILDLIBS) $(GPTFLAGS) -o $@

$(IOBJS): $(SRCS)
	$(ICC) $(IFLAGS) $^ -c

%_g.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

cleanall:
	$(RM) $(OBJS) $(IOBJS) $(TARGETS) *~

