TARGETS = release debug
CC = gcc
CFLAGS = -Wall -Wextra
LDLIBS= -lm -lopenblas
SRCS = my_dense.c my_sparse.c timer.c csr.c spmv.c
OBJS = $(SRCS:.c=.o)
OBJS_DBG = $(SRCS:.c=_dbg.o)

all: $(TARGETS)

%.o: %.c
	$(CC) $(CFLAGS) -O2 -c $<

%_dbg.o: %.c
	$(CC) $(CFLAGS) -g -O0 -c -o $@ $<

release: $(OBJS)
debug: $(OBJS_DBG)

$(TARGETS):
	$(CC) $^ $(LDFLAGS) $(LDLIBS) -o $@
