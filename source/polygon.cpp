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
    int i;
    for(i = 0; i < shape.size() - 1; i++) {
        this->points.push_back(points[shape[i]]);
        int *tmp = new int [2];
        tmp[0] = shape[i];
        tmp[1] = shape[i + 1];
        segments.push_back(tmp);
    }
    this->points.push_back(points[shape[i]]);
    int *tmp = new int [2];
    tmp[0] = shape[i];
    tmp[1] = shape[0];
    segments.push_back(tmp);
    this->shape = shape;
    find_perimeter();
    find_id();
}

/* constructor: initialize with segments */
Polygon::Polygon(std::vector<int> shape, std::vector<int *> segments, Point *points) {
    int i;
    for(i = 0; i < segments.size(); i++) {
        bool found = false;
        for(Point p : this->points) {
            if(p.index == points[segments[i][0]].index) {
                found = true;
                break;
            }
        }
        if(!found) {
            this->points.push_back(points[segments[i][0]]);
        }
        found = false;
        for(Point p : this->points) {
            if(p.index == points[segments[i][1]].index) {
                found = true;
                break;
            }
        }
        if(!found) {
            this->points.push_back(points[segments[i][1]]);
        }
    }
    this->shape = shape;
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

/* tests if a point is part of the polygon
 * @param P, the point to test
 * @return the index of the point if found, -1 otherwise
 */
int Polygon::point_match(Point P) {
    for(int i = 0; i < points.size(); i++) {
        if(P.equals(points[i])) {
            return i;
        }
    }

    return -1;
}

/* tests if a segment is part of the polygon
 * @param beginning, the start index of the segment to test
 * @param end, the end index of the segment to test
 * @return the index of the segment if found, -1 otherwise
 */
int Polygon::segment_match(int beginning, int end) {
    for(int i = 0; i < segments.size(); i++) {
        if(((segments[i][0] == beginning) && (segments[i][1] == end)) || ((segments[i][0] == end) && (segments[i][1] == beginning))) {
            return i;
        }
    }

    return -1;
}
