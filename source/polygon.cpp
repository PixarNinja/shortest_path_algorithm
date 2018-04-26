/*
 * Implementation of the Vector class and methods used for vector analysis
 *
 * Shortest Path Algorithm
 * Mark Wesley Harris
 * 2018
 */

#include <vector>
#include <new>
#include "polygon.h"

Polygon::Polygon() {}

/* constructor: base */
Polygon::Polygon(std::vector<int> shape, Point *points) {
    this->shape = shape;
    int i;
    for(i = 0; i < shape.size() - 1; i++) {
        this->points.push_back(points[shape[i]]);
        int *tmp = new int [2];
        tmp[0] = shape[i];
        tmp[1] = shape[i + 1];
        segments.push_back(tmp);
    }
    int *tmp = new int [2];
    tmp[0] = shape[i];
    tmp[1] = shape[0];
    segments.push_back(tmp);
    find_perimeter();
    find_id();
}

/* constructor: clone a Polygon */
Polygon::Polygon(const Polygon &S) {
    shape = S.shape;
    segments = S.segments;
    points = S.points;
    perimeter = S.perimeter;
    id = S.id;
}
