CC=gcc
CPP=g++
CFLAGS=-g
CPPFLAGS=-g -std=c++11

all: brute_force compare plot shortest_path

brute_force: brute_force.c
	$(CC) $(CFLAGS) brute_force.c -o brute_force -lm

compare: compare.c
	$(CC) $(CFLAGS) compare.c -o compare

multiple: multiple_contours.cpp
	$(CPP) $(CPPFLAGS) multiple_contours.cpp -o multiple_contours -lm

plot: plot.c
	$(CC) $(CFLAGS) plot.c -o plot

pop: pop_reverse.cpp
	$(CPP) $(CPPFLAGS) pop_reverse.cpp -o pop_reverse -lm

tao2: tao2.c construct_contour.c
	$(CC) $(CPPFLAGS) construct_contour.c -c -lm
	$(CC) $(CPPFLAGS) tao2.c construct_contour.o -o tao2 -lm

tessellate: tessellate.cpp
	$(CPP) $(CPPFLAGS) tessellate.cpp -o tessellate -lm

shortest_path: shortest_path.cpp
	$(CPP) $(CPPFLAGS) shortest_path.cpp -o shortest_path -lm

clean: 
	rm brute_force compare multiple_contours plot pop_reverse tao2 tessellate
