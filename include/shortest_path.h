/*
 * Definition of shortest_path.cpp functions for main()
 *
 * Shortest Path Algorithm
 * Mark Wesley Harris
 * 2018
 */

#ifndef SHORTESTPATH_H
#define SHORTESTPATH_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <ctype.h>
#include <unistd.h>

#include <new>
#include <vector>
#include <iostream>
#include <deque>
#include <algorithm>
#include <iterator>
#include <array>
#include <string>
#include <unordered_map>
#include <sstream>

#include "polygon.h"
#include "point.h"
#include "vector.h"

using namespace std;

/* print functions */
void memory_error(void);

/* index matching functions */
int edge_match(Polygon polygon, int *edge);
int segment_match(vector<int *> segments, int beginning, int end);
int shape_match(vector<int> shape, int vertex);
int point_match(Point *points, int size, int vertex);

/* segment matching functions */
vector<int *> edge_search(vector<int *> segments, int vertex, Point *points, int size);
vector<int *> disjoint_edges(Polygon A, Polygon B);
vector<int *> shared_edges(Polygon A, Polygon B);

/* point matching functions */
vector<int> shared_points(Polygon A, Polygon B);

/* edge search functions */
vector<int> breadth_first_index_search(vector<int *> segments, vector<int> *processed, vector<int> seed, Point *points, int size);
vector<int> index_search(vector<int *> segments, vector<int> *processed, int vertex, Point *points, int size);

/* boolean check functions */
bool intersection(Vector V1, Vector V2);
bool same_direction(Vector V1, Vector V2);

/* polygon processing functions */
Polygon add_polygons(Polygon A, Polygon B, Point *points);
Polygon sub_polygons(Polygon A, Polygon B, Point *points);

/* math calculation functions */
double angle(Vector V1, Vector V2);
double distance_p(Point P1, Point P2);
double distance_v(Vector V1, Vector V2);
double dot_product(Vector V1, Vector V2);
Vector projection(Vector V1, Vector V2);
Vector circular_gradient(Point center, double radius, Point p);
Point find_intersection(Vector V1, Vector V2);
double determinant(Vector V1, Vector V2);
bool overlap(Vector V1, Vector V2);
Polygon find_convex_hull(Point *points, int size);

/* shortest path algorithm functions */
vector<int *> all_w_segments(Point *points, int size);
bool test_w_segment(vector<int *> segments, Vector L, double interval, Point *points, int n);
vector<Polygon> generate_final_paths(vector<int *> segments, Point *points, int n);
Polygon dijkstra_polygon(vector<int *> segments, int *segment, Point *points, int n);

#endif
