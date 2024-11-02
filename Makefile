#module load cesga/2022 iimkl
TARGETS = gsl mkl
CC = gcc
ICC = icc
OFLAGS = -O0
CFLAGS = -Wall -Wextra -D_GSL_
IFLAGS = -g3 -Wall -Wextra -D_MKL_ -I"${MKLROOT}/include" -diag-disable=10441
LDLIBS = -lm -lgsl -lopenblas
ILDLIBS = -lm -qmkl
#ILDLIBS = -L"${MKLROOT}/lib/intel64" -lm -qmkl

ISRCS = timer.c my_dense.c my_csr.c my_csc.c my_coo.c spmv.c
SRCS = timer.c my_dense.c my_csr.c my_coo.c my_csc.c spmv.c
IOBJS = $(ISRCS:.c=.o)
OBJS = $(SRCS:.c=_g.o)

all: $(TARGETS)

gsl: $(OBJS) 	
	$(CC) $(OFLAGS) $(CFLAGS) $^ $(LDLIBS) -o $@
mkl: $(IOBJS)
	$(ICC) $(OFLAGS) $(IFLAGS) $^ $(ILDLIBS) -o $@

$(IOBJS): $(ISRCS)
	$(ICC) $(IFLAGS) $^ -c

%_g.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

cleanall:
	$(RM) $(OBJS) $(IOBJS) release *~

#cleanall:
#	$(RM) $(OBJS) $(OBJS_DBG) debug release *~

