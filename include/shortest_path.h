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

#include "point.h"
#include "vector.h"

#define NUM_FILES 8

using namespace std;

struct polygon_t {
    vector<int> shape;
    double perimeter;
};

void construct_segments(vector<int *> *segments, Point *points, Point begin, int n, int size, FILE *gnu_files[NUM_FILES], int *mapped, int **recorded);
void join_vertex(vector<int *> *segments, Point *points, Point begin, int n, int size);
void join_segment(vector<int *> *segments, Point *points, Point begin, Point end, int n, int m, int size);
vector<struct polygon_t> construct_polygons(vector<int *> segments, Point *points, int size);
vector<vector<int> > tessellate(vector<vector<int> > tessellations, vector<int *> segments, Point *points, int size, char init, char add, int branch);
vector<struct polygon_t> tessellate_cross(vector<int *> segments, int i, int j, Point *points, int size);
vector<int> *find_shape(vector<int *> segments, Point *points, Point start, int size, char init, char add, Vector X, Vector Y);
vector<int> init_path(vector<int> path, vector<int *> edges, Point *points, Vector X, Vector Y, char type, int branch);
vector<int> add_path(vector<int> path, vector<int *> edges, Point *points, Vector X, Vector Y, char type);
vector<int *> edge_search(vector<int *> segments, int vertex, Point *points, int size);
int index_match(vector<int *> segments, int vertex);
int shape_search(vector<int> shape, int vertex);
int edge_match(struct polygon_t polygon, int *edge);
int polygons_search(vector<vector<int> > polygons, int vertex);
double find_perimeter(vector<int> shape, Point *points);
int segment_match(vector<int *> segments, int beginning, int end);
int duplicate_search(vector<int> shape);
vector<struct polygon_t> delete_duplicates(vector<struct polygon_t> polygons);
vector<struct polygon_t> optimize_polygons(vector<struct polygon_t> polygons, vector<int *> *segments, Point *points, int size);
void remove_crosses(vector<int *> *segments, Point *points, int size);
void finalize_segments(vector<int *> *segments, Point *points, int size);
bool intersection(Vector V1, Vector V2);
struct polygon_t find_shortest_path(vector<struct polygon_t> polygons, Point *points, int size);
int accept_polygon(struct polygon_t polygon, vector<int *> segments, Point *points);
int smallest_neighbour(vector<struct polygon_t> polygons, struct polygon_t source, int n);
vector<int *> disjoint_edges(struct polygon_t A, struct polygon_t B);
vector<int *> shared_edges(struct polygon_t A, struct polygon_t B);
vector<int> shared_points(struct polygon_t A, struct polygon_t B);
void visit_polygon(int *visited, struct polygon_t polygon, Point *points);
struct polygon_t add_polygons(struct polygon_t A, struct polygon_t B, Point *points);
struct polygon_t sub_polygons(struct polygon_t A, struct polygon_t B, Point *points);
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
vector<int *> midpoint_construction(Point *points, int size, FILE *gnu_files[NUM_FILES]);
bool test_w_segment(Vector L, double interval, Point *points, int n);
double determinant(Vector V1, Vector V2);
vector<Point> generate_w_points(vector<Point> w_points, Vector L, double interval);
Point minimum_tao_distance(Vector V, Point *points, int size);
vector<int *> remove_crossing_segments(vector<int *> segments, int s, Point *points);
vector<int *> fix_overlap(int *test, vector<int *> segments, Point *points);
struct polygon_t create_polygon(int *edge, vector<int *> segments, Point *points, int size);

#endif
