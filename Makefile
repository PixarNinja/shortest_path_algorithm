CC=gcc
CPP=g++
CFLAGS=-g

all: multiple shortest tao2 pop tessellate

shortest: shortest.c
	$(CC) $(CFLAGS) shortest.c -o shortest -lm

tao2: tao2.c construct_contour.c
	$(CC) $(CFLAGS) construct_contour.c -c -lm
	$(CC) $(CFLAGS) tao2.c construct_contour.o -o tao2 -lm

multiple: multiple_contours.cpp
	$(CPP) $(CFLAGS) multiple_contours.cpp -o multiple_contours -lm

pop: pop_reverse.cpp
	$(CPP) $(CFLAGS) pop_reverse.cpp -o pop_reverse -lm

tessellate: tessellate.cpp
	$(CPP) $(CFLAGS) tessellate.cpp -o tessellate -lm

clean: 
	rm shortest multiple_contours tao2 pop_reverse tessellate
