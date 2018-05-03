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

#include "polygon.h"
#include "point.h"
#include "vector.h"

#define NUM_FILES 8

using namespace std;

vector<int *> edge_search(vector<int *> segments, int vertex, Point *points, int size);
int index_match(vector<int *> segments, int vertex);
int shape_search(vector<int> shape, int vertex);
int edge_match(Polygon polygon, int *edge);
int polygons_search(vector<vector<int> > polygons, int vertex);
double find_perimeter(vector<int> shape, Point *points);
int segment_match(vector<int *> segments, int beginning, int end);
int duplicate_search(vector<int> shape);
bool intersection(Vector V1, Vector V2);
bool same_direction(Vector V1, Vector V2);
Polygon find_shortest_path(vector<Polygon> polygons, Point *points, int size);
int smallest_neighbour(vector<Polygon> polygons, Polygon source, int n);
vector<int *> disjoint_edges(Polygon A, Polygon B);
vector<int *> shared_edges(Polygon A, Polygon B);
vector<int> shared_points(Polygon A, Polygon B);
void visit_polygon(int *visited, Polygon polygon, Point *points);
Polygon add_polygons(Polygon A, Polygon B, Point *points);
Polygon sub_polygons(Polygon A, Polygon B, Point *points);
int point_match(Point *points, int size, int vertex);
double calculate_curvature(Vector T1, Vector T2, double tao);
double angle(Vector V1, Vector V2);
double angle_t(double tao);
double tao_distance(Vector V, double curvature, double theta);
double distance_p(Point P1, Point P2);
double distance_v(Vector V1, Vector V2);
double dot_product(Vector V1, Vector V2);
void print(Vector V, Vector T1, Vector T2, double curvature, double theta, double tao, double tao_distance);
void memory_error(void);
vector<Polygon> construct_w_polygons(Polygon base, Point *points, int size, FILE *gnu_files[NUM_FILES]);
vector<Polygon> tessellate_w_polygon(Polygon S, double interval, Point *points, int size);
bool test_w_segment(Vector L, double interval, Point *points, int n);
double determinant(Vector V1, Vector V2);
vector<Point> generate_w_points(vector<Point> w_points, Vector L, double interval);
Point minimum_tao_distance(Vector V, Point *points, int size);
vector<int *> remove_crossing_segments(vector<int *> segments, int s, Point *points, int size);
vector<int *> fix_overlap(int *test, vector<int *> segments, Point *points);
bool overlap(Vector V1, Vector V2);
vector<Polygon> create_polygon(int *edge, vector<int *> segments, Point *points, int size);
vector<Polygon> delete_duplicate_polygons(vector<Polygon> polygons, Point *points);
Polygon find_convex_hull(Point *points, int size);

#endif
