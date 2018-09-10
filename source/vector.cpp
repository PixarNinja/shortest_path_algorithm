/*
 * Implementation of the Vector class and methods used for vector analysis
 *
 * Shortest Path Algorithm
 * Mark Wesley Harris
 * 2018
 */

#include "vector.h"
#include "point.h"
#include <string>
#include <iostream>

Vector::Vector() {}

Vector::Vector(std::string name, Point start, Point end) {
    this->name = name;
    this->start = start;
    this->end = end;
    refresh();
}

/* creates a normalized vector centered at start
 * and facing end
 */
Vector::Vector(std::string name, Point start, Point end, int index) {
    double x_offset = end.x - start.x;
    double y_offset = end.y - start.y;
    this->name = name;
    this->start = Point(end);
    /* create a normalized point centered at start in the direction of end */
    this->end = Point(this->start.x + (x_offset / distance_p(start, end)), this->start.y + (y_offset / distance_p(start, end)), index);
    refresh();
}

/* normalizes the vector */
void Vector::normalize() {
    end = Point(start.x + (i / length), start.y + (j / length), -1);
    refresh();
}

/* refreshes i and j, and stores the length of the vector */
void Vector::refresh() {
    i = end.x - start.x;
    j = end.y - start.y;
    length = sqrt(pow(i, 2) + pow(j, 2));
}

/* scoots the vector to start at its endpoint */
void Vector::scoot() {
    /* shift points */
    start = Point(end);
    end.offset(i, j);
}

/* offsets a vector
 * @param x_offset, the x offset
 * @param y_offset, the y offset
 */
void Vector::offset(double x_offset, double y_offset) {
    start.offset(x_offset, y_offset);
    end.offset(x_offset, y_offset);
}

/* tests if a vector is equal
 * @param V, the vector to test
 * @return true if equal, false otherwise
 */
bool Vector::equals(Vector V) {
    if((start.equals(V.start) && end.equals(V.end)) || (start.equals(V.end) && end.equals(V.start))) {
        return true;
    }
    return false;
}

/* prints information about this vector */
void Vector::print() {
    std::cout << name;
    printf(": datapoints = %d (%0.3lf,%0.3lf), %d (%0.3lf,%0.3lf)\n   components = <%0.3lf,%0.3lf>\n   length     = %0.3lf\n\n", start.index, start.x, start.y, end.index, end.x, end.y, i, j, length);
}
