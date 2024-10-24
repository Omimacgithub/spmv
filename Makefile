TARGETS = release debug
CC = gcc
CFLAGS = -Wall -Wextra
LDLIBS= -lm -lgsl -lopenblas

SRCS = my_dense.c my_csr.c my_csc.c my_coo.c timer.c spmv.c
OBJS = $(SRCS:.c=.o)
OBJS_DBG = $(SRCS:.c=_dbg.o)

all: $(TARGETS)

%.o: %.c
	$(CC) $(CFLAGS) -O2 -c $<

%_dbg.o: %.c
	$(CC) $(CFLAGS) -O0 -g3 -c -o $@ $<

release: $(OBJS)
debug: $(OBJS_DBG)

$(TARGETS):
	$(CC) $^ $(LDFLAGS) $(LDLIBS) -o $@

clean:
	$(RM) $(OBJS) $(OBJS_DBG) *~

cleanall:
	$(RM) $(OBJS) $(OBJS_DBG) debug release *~

