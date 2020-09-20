CC=gcc

CFLAGS= -I./src/include

AR = ar
RANLIB = ranlib

SRCS = src/mklist.c 
OBJS = $(SRCS:.c=.o)

all	: mklist

mklist	: $(OBJS)
	$(CC) $(CFLAGS) -o mklist $(OBJS)


