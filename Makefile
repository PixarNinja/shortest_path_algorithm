CC=gcc
CPP=g++
CFLAGS=-g
CPPFLAGS=-g -std=c++11

all: tessellate

shortest: shortest.c
	$(CC) $(CPPFLAGS) shortest.c -o shortest -lm

tao2: tao2.c construct_contour.c
	$(CC) $(CPPFLAGS) construct_contour.c -c -lm
	$(CC) $(CPPFLAGS) tao2.c construct_contour.o -o tao2 -lm

multiple: multiple_contours.cpp
	$(CPP) $(CPPFLAGS) multiple_contours.cpp -o multiple_contours -lm

pop: pop_reverse.cpp
	$(CPP) $(CPPFLAGS) pop_reverse.cpp -o pop_reverse -lm

tessellate: tessellate.cpp
	$(CPP) $(CPPFLAGS) tessellate.cpp -o tessellate -lm

clean: 
	rm shortest multiple_contours tao2 pop_reverse tessellate
