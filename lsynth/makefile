CC=gcc

CFLAGS=-ggdb

AR = ar
RANLIB = ranlib

## This is lame.  I don't know how to do both .c and .cpp for the OBJS line
## so I pretend L3*.cpp is L3*.c and make rules for them later.
#
SRCS = src/lsynthcp.c src/band.c src/hose.c src/curve.c src/mathlib.c
OBJS = $(SRCS:.c=.o)

all	: lsynthcp

lsynthcp:   $(OBJS) $(LIBS)
	$(CC) $(CFLAGS) $(OBJS) -lm -o lsynthcp
