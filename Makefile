PROG = test

OBJS = test.o corpsfinis.o

CC = gcc
CPFLAGS = -Wall -Wpointer-arith -O2


test: test.o corpsfinis.o
	$(CC) $(CPFLAGS) -o test  test.o corpsfinis.o

test.o: corpsfinis.h test.h test.c
	$(CC) $(CPFLAGS) -c test.c

corpsfinis.o: corpsfinis.c corpsfinis.h
	$(CC) $(CPFLAGS) -c corpsfinis.c

clean:
	rm -f $(OBJS)


