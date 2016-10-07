CC=gcc
CFLAGS=-g -Wall

all: tao2 shortest

shortest:
	$(CC) $(CFLAGS) shortest.c -o shortest -lm

tao2:
	$(CC) $(CFLAGS) tao2.c -o tao2 -lm

clean: 
	rm shortest tao2
