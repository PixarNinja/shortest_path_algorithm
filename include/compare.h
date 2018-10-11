/*
 * Definition of compare.cpp functions for main()
 *
 * Shortest Path Algorithm
 * Mark Wesley Harris
 * 2018
 */

#ifndef COMPARE_H
#define COMPARE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string>
#include <vector>

#include "point.h"

double distance(Point P1, Point P2);
int point_match(Point *points, int size, int vertex);

using namespace std;

#endif
