CC=gcc
CPP=g++
CFLAGS=-g -Wall

all: multiple shortest tao2

shortest:
	$(CC) $(CFLAGS) shortest.c -o shortest -lm

tao2:
	$(CC) $(CFLAGS) construct_contour.c -c -lm
	$(CC) $(CFLAGS) tao2.c construct_contour.o -o tao2 -lm

multiple:
	$(CPP) $(CFLAGS) multiple_contours.cpp -o multiple_contours -lm

clean: 
	rm shortest multiple tao2
