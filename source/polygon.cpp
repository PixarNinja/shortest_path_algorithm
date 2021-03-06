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
    this->segments = segments;
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

    /* create remaining point vector */
    std::vector<Point> remaining;
    for(i = 0; i < points.size(); i++) {
        if(i == k) {
            continue;
        }
        remaining.push_back(Point(points[i]));
    }
    remaining.push_back(Point(points[k])); // ensure points[k] is not first

    /* create initial segment by searching for the least polar angle */
    Vector M = Vector("M", points[k], remaining[0]);
    int m = 0;
    for(i = 1; i < remaining.size() - 1; i++) {
        Vector X = Vector("X", points[k], points[k]);
        X.end.offset(1.0, 0.0);
        X.refresh();
        Vector V = Vector("V", points[k], remaining[i]);

        /* find the minimum polar angle from points[k] */
        double curr = angle(X, V) * 180 / M_PI; // we know that det(X, V) >= 0
        double min = angle(X, M) * 180 / M_PI; // we know that det(X, M) >= 0
        if(curr < min) {
            M = Vector(V);
            m = i;
        }
        else if(curr == min) { // if equal check the vector length
            if(V.length < M.length) {
                M = Vector(V);
                m = i;
            }
        }
    }
    
    /* remove M.end from remaining vector */
    remaining.erase(remaining.begin() + m);

    /* push starting point */
    std::vector<int> shape;
    shape.push_back(point_match(data, size, points[k].index));
    shape.push_back(point_match(data, size, M.end.index));

    /* loop until we reach points[k] */
    while(remaining.size() > 0 && !M.end.equals(points[k])) {
        /* normalize and shift */
        Vector T1 = Vector(M);
        T1.start = T1.end;
        T1.end.offset(T1.i, T1.j);
        T1.normalize();
        T1.refresh();

        /* find the next segment */
        M = Vector("M", T1.start, remaining[0]);
        m = 0;
        for(i = 1; i < remaining.size(); i++) {
            /* find the minimum polar angle from points[k] */
            Vector T2 = Vector("T2", T1.start, remaining[i]);
            double curr = angle(T1, T2) * 180 / M_PI; // we know that det(T1, T2) >= 0
            double min = angle(T1, M) * 180 / M_PI; // we know that det(T1, M) >= 0
            if(curr < min) {
                M = Vector(T2);
                m = i;
            }
            else if(curr == min) { // if equal check the vector length
                if(T2.length < M.length) {
                    M = Vector(T2);
                    m = i;
                }
            }
        }

        /* remove M.end from remaining vector */
        remaining.erase(remaining.begin() + m);

        /* push found point */
        shape.push_back(point_match(data, size, M.end.index));
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
    bool run_test = false;

    /* first pass, see if the point should be tested */
    for(i = 1; i < shape.size(); i++) {
        int start = shape[i - 1];
        int end = shape[i];
        Vector V = Vector("V", points[start], points[end]);
        Vector T = Vector("T", points[start], test);
        Vector P = projection(T, V);
        if(P.length <= V.length) {
            run_test = true;
            break;
        }
    }
    /* create a segment from each point in the shape */
    if(run_test) {
        for(i = 1; i < shape.size(); i++) {
            int start = shape[i - 1];
            int end = shape[i];
            Vector V = Vector("V", points[start], points[end]);
            Vector T = Vector("T", points[start], test);
            Vector P = projection(T, V);
            if(P.length <= V.length) {
                /* test if the point is on the left side or in-line */
                if(determinant(V, T) < 0) {
                    return false;
                }
            }
        }

        return true;
    }

    return false;
}

/* tests if a polygon is equal
 * @param S, the polygon to test
 * @return true if equal, false otherwise
 */
bool Polygon::equals(Polygon S) {
    if(id == S.id) {
        return true;
    }
    return false;
}
