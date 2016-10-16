CC=gcc
CFLAGS=-g -Wall

all: tao2 shortest

shortest:
	$(CC) $(CFLAGS) shortest.c -o shortest -lm

tao2:
	$(CC) $(CFLAGS) construct_contour.c -c -lm
	$(CC) $(CFLAGS) tao2.c construct_contour.o -o tao2 -lm

clean: 
	rm shortest tao2
