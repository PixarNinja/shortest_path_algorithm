#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <ctype.h>
#include <unistd.h>

struct point_t {
    double x;
    double y;
    double tao_d;
    int index;
};

struct vector_t {
    struct point_t point[2];
    double length;
    double i;
    double j;
};

struct curvature_t {
    struct vector_t T1;
    struct vector_t T2;
    struct vector_t V;
    double tao;
    double theta;
    double curvature;
};


int **construct_contour(struct point_t *points, int size);
double calculate_curvature(struct curvature_t k);
double calculate_theta(struct curvature_t k);
double tao_distance(struct curvature_t k);
double angle_v(struct vector_t start, struct vector_t end);
double distance_p(struct point_t start, struct point_t end);
double distance_v(struct vector_t start, struct vector_t end);
double length_v(struct vector_t v);
double dot_product(struct vector_t start, struct vector_t end);
void print_k(struct curvature_t k);
