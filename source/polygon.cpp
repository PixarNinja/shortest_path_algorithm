/*
 * Implementation of the Vector class and methods used for vector analysis
 *
 * Shortest Path Algorithm
 * Mark Wesley Harris
 * 2018
 */

#include <vector>
#include "polygon.h"

Polygon::Polygon() {}

/* constructor: base */
Polygon::Polygon(std::vector<int> shape, Point *points) {
    this->shape = shape;
    for(int index : shape) {
        this->points.push_back(points[index]);
    }
    find_perimeter();
    find_id();
}

/* constructor: clone a Polygon */
Polygon::Polygon(const Polygon &S) {
    shape = S.shape;
    points = S.points;
    perimeter = S.perimeter;
    id = S.id;
}
