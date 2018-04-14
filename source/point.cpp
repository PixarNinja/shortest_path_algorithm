/*
 * Implementation of the Vector class and methods used for vector analysis
 *
 * Shortest Path Algorithm
 * Mark Wesley Harris
 * 2018
 */

#include "vector.h"
#include <string>

Point::Point() {}

/* constructor: base */
Point::Point(double x, double y, int index) {
    this->x = x;
    this->y = y;
    this->index = index;
    theta = 0;
    curvature = 0;
    tao = 0;
    tao_distance = 0;
}

/* constructor: verbose */
Point::Point(double x, double y, int index, double theta, double curvature, double tao, double tao_distance) {
    this->x = x;
    this->y = y;
    this->index = index;
    this->theta = theta;
    this->curvature = curvature;
    this->tao = tao;
    this->tao_distance = tao_distance;
}


/* constructor: clone a Point */
Point::Point(const Point &P) {
    x = P.x;
    y = P.y;
    theta = P.theta;
    curvature = P.curvature;
    tao = P.tao;
    tao_distance = P.tao_distance;
    index = P.index;
}

/* offsets the point
 * @param x, the x offset
 * @param y, the y offset
 */
void Point::offset(double x_offset, double y_offset) {
    x += x_offset;
    y += y_offset;
}

/* tests if a point is equal
 * @param P, the point to test
 * @return true if equal, false otherwise
 */
bool Point::equals(Point P) {
    if(index == P.index) {
        return true;
    }
    return false;
}
