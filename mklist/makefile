CC=gcc

CFLAGS= -I./include

AR = ar
RANLIB = ranlib

SRCS = mklist.c 
OBJS = $(SRCS:.c=.o)

all	: mklist

mklist	: $(OBJS)
	$(CC) $(CFLAGS) -o mklist $(OBJS)


