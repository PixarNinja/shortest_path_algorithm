/*
 * Definition of brute_force.cpp functions for main()
 *
 * Shortest Path Algorithm
 * Mark Wesley Harris
 * 2018
 */

#ifndef BRUTEFORCE_H
#define BRUTEFORCE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <vector.h>
#include <algorithm>

#include "point.h"
#include "vector.h"

using namespace std;

#define SIZE 7

struct data_t {
    Vector v;
    Point p;
    int i;
    double curr;
} point;

double angle(Vector V1, Vector V2);
double dot_product(Vector V1, Vector V2);
bool intersection(Vector V1, Vector V2);
bool same_direction(Vector V1, Vector V2);
double determinant(Vector V1, Vector V2);
void swap(struct data_t *x, struct data_t *y);
void permute_dist(struct data_t *point, struct data_t *shortest, int index, int n);
void permute_curv(struct data_t *point, struct data_t *shortest, int index, int n);
void permute_tao(struct data_t *point, struct data_t *shortest, int index, int n);
double distance(struct data_t *point, int first, int last);
double curvature(struct data_t *point, int prev, int curr, int next, vector<int *> path);
double tao(struct data_t *point, int prev, int curr, int next, vector<int *> path);
double distance_p(Point P1, Point P2);
int factorial(int n);

#endif
