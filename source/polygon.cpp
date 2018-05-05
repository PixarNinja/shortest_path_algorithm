/*
 * Implementation of the Vector class and methods used for std::vector analysis
 *
 * Shortest Path Algorithm
 * Mark Wesley Harris
 * 2018
 */

#include <vector>
#include <new>
#include <cfloat>
#include "polygon.h"
#include "vector.h"

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

/* tests if a point is part of a point array
 * @param points, the array of points to test
 * @param size, the size of the point array
 * @param vertex, the vertex value to find
 * @return the index of the point if found, -1 otherwise
 */
int Polygon::point_match(Point *points, int size, int vertex) {
    for(int i = 0; i < size; i++) {
       if(points[i].index == vertex) {
            return i;
        }
    }
    return -1;
}

/* tests if a segment is part of a segment vector
 * @param segments, the segement vector to test
 * @param beginning, the start index of the segment to test
 * @param end, the end index of the segment to test
 * @return the index of the segment if found, -1 otherwise
 */
int Polygon::segment_match(std::vector<int *> segments, int beginning, int end) {
    for(int i = 0; i < segments.size(); i++) {
        if(((segments[i][0] == beginning) && (segments[i][1] == end)) || ((segments[i][0] == end) && (segments[i][1] == beginning))) {
            return i;
        }
    }

    return -1;
}

/* re-orders the points of the shape into a hull
 * @param data, the overall point array
 * @param size, the size of the array
 */
void Polygon::create_hull(Point *data, int size) {
    int i;
    int j;

    /* find the lowest, right-most point */
    double min_y = DBL_MAX;
    double max_x = DBL_MIN;
    int k = -1;
    for(i = 0; i < points.size(); i++) {
        if(points[i].y <= min_y) {
            /* if the y value is equal to min, check x value */
            if(points[i].y == min_y) {
                if(points[i].x > max_x) {
                    k = i;
                    max_x = points[i].x;
                }
            }
            /* else it is less than the min, so reset x value */
            else {
                k = i;
                max_x = DBL_MIN;
            }
            min_y = points[i].y;
        }
    }

    /* create remaining point std::vector */
    std::vector<Point> remaining;
    for(i = 0; i < points.size(); i++) {
        if(i != k) {
            remaining.push_back(Point(points[i]));
        }
    }

    /* bubble sort the remaining points by polar-angle */
    Vector X = Vector("X", points[k], points[k]);
    X.end.offset(1.0, 0.0);
    X.refresh();
    for(i = 0; i < remaining.size(); i++) {
        for(j = remaining.size() - 1; j > i; j--) {
            Vector A = Vector("A", points[k], remaining[j]);
            Vector B = Vector("B", points[k], remaining[j - 1]);
            /* find the polar angle */
            double curr;
            if(determinant(X, A) >= 0) {
                curr = angle(A, X) * 180 / M_PI;
            }
            else {
                curr = 360 - (angle(A, X) * 180 / M_PI);
            }
            double prev;
            if(determinant(X, B) >= 0) {
                prev = angle(B, X) * 180 / M_PI;
            }
            else {
                prev = 360 - (angle(B, X) * 180 / M_PI);
            }
            if(curr < prev) {
                Point tmp = Point(remaining[j]);
                remaining[j] = remaining[j - 1];
                remaining[j - 1] = tmp;
            }
            /* if equal check the std::vector length */
            else if(curr == prev) {
                /* TODO: FIX! check if they lie on a horizontal line */
                if(A.j == 0.0 && B.j == 0.0) {
                    if(A.length > B.length) {
                        Point tmp = Point(remaining[j]);
                        remaining[j] = remaining[j - 1];
                        remaining[j - 1] = tmp;
                    }
                }
                else {
                    if(A.length < B.length) {
                        Point tmp = Point(remaining[j]);
                        remaining[j] = remaining[j - 1];
                        remaining[j - 1] = tmp;
                    }
                }
            }
        }
    }

    /* push starting values onto the stack */
    std::vector<Point> stack;
    stack.push_back(points[k]);
    stack.push_back(remaining[0]);
    stack.push_back(remaining[1]);

    /* search for non-left turns and pop them off */
    for(i = 2; i < remaining.size(); i++) {
        Vector V = Vector("V", stack[stack.size() - 1], stack[stack.size() - 2]); // current line segment
        Vector T = Vector("T", stack[stack.size() - 1], remaining[i]); // test line segment
        while(determinant(V, T) > 0) { // look for left turns
            printf("POPPING: %d\n", stack.back().index);
            stack.erase(stack.begin() + stack.size() - 1);
            if(stack.size() <= 2) {
                break;
            }
            V = Vector("V", stack[stack.size() - 1], stack[stack.size() - 2]); // current line segment
            T = Vector("T", stack[stack.size() - 1], remaining[i]); // test line segment
        }
        stack.push_back(remaining[i]);
    }
    stack.push_back(points[k]);

    /* create shape off of original points (function parameter) */
    std::vector<int> shape;
    for(Point p : stack) {
        shape.push_back(point_match(data, size, p.index));
    }

    /* store new shape */
    this->shape = shape;
}

/* finds if a point is inside of the shape (requires the shape to be sorted as a hull)
 * @param test, the point to test
 * @param points, the point array in which the test point exists
 * @return true if the point is inside the shape, false otherwise
 */
bool Polygon::contains(Point test, Point *points, int size) {
    int i = 0;

    /* create a segment from each point in the shape */
    for(i = 1; i < shape.size(); i++) {
        Vector V = Vector("V", points[shape[i - 1]], points[shape[i]]);
        Vector T = Vector("T", points[shape[i - 1]], test);

        /* test if the point is on the left side or in-line */
        if(determinant(V, T) < 0) {
            return false;
        }
    }

    return true;
}
