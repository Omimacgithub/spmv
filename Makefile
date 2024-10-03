TARGETS = spmv
CC = gcc
CFLAGS = -Wall -Wextra
LDLIBS= -lm -lopenblas

SRCS = my_dense.c my_sparse.c timer.c spmv.c
OBJS = $(SRCS:.c=.o)

all: $(TARGETS)

%.o: %.c
	$(CC) $(CFLAGS) -O2 -c $<

spmv: $(OBJS)

$(TARGETS):
	$(CC) $^ $(LDFLAGS) $(LDLIBS) -o $@
