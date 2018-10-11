/*
 * Shortest path algorithm function definitions
 *
 * Shortest Path Algorithm
 * Mark Wesley Harris
 * 2018
 */

#include "shortest_path.h"

#ifndef SHORTESTPATH_H

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

#endif

/* prints to the terminal if there is an error assigning memory */
void memory_error(void)
{
    printf("\n\nError assigning memory. Exiting Program. Good Day.\n\n");
}

/* searches through a shape for a matching edge
 * @param polygon, the shape to search in
 * @param edge, the edge to search for
 * @return the index of the edge or -1 if not found
 */
int edge_match(Polygon polygon, int *edge)
{
    vector<int *> segments;
    int *tmp;
    for(int i = 0; i < polygon.shape.size() - 1; i++) {
        tmp = new int [2];
        tmp[0] = polygon.shape[i];
        tmp[1] = polygon.shape[i + 1];
        segments.push_back(tmp);
    }
    return segment_match(segments, edge[0], edge[1]);
}

/* searches through a vector of segments for a matching segment
 * @param segments, all edges which form graph G
 * @param beginning, the beginning index
 * @param ending, the ending index
 * @return the index of the edge or -1 if not found
 */
int segment_match(vector<int *> segments, int beginning, int end)
{
    /* check if the segment is found */
    for(int i = 0; i < segments.size(); i++) {
        if(((segments[i][0] == beginning) && (segments[i][1] == end)) || ((segments[i][0] == end) && (segments[i][1] == beginning))) {
            return i;
        }
    }
    return -1;
}

/* searches through a shape for a matching vertex
 * @param shape, the shape to search through
 * @param vertex, the vertex to search for
 * @return the index of the edge or -1 if not found
 */
int shape_match(vector<int> shape, int vertex)
{
    for(int i = 0; i < shape.size() - 1; i++) {
        /* check if the vertex is found */
        if(shape[i] == vertex) {
            return i;
        }
    }
    return -1;
}

/* returns the index of the requested vertex */
int point_match(Point *points, int size, int vertex)
{
    for(int i = 0; i < size; i++) {
       if(points[i].index == vertex) {
            return i;
        }
    }
    return -1;
}

/* searches through a vector of segments for all matching segment
 * @param segments, all edges which form graph G
 * @param vertex, the Point.index value to search for
 * @param points, the dataset of points
 * @param size, the size of the dataset
 * @return edges, a vector of matching segment elements
 */
vector<int *> edge_search(vector<int *> segments, int vertex, Point *points, int size)
{
    vector<int *> edges;
    int *tmp;
    int i = 0;
    int j = 0;
    int index;

    /* find the indexed value */
    for(index = 0; index < size; index++) {
        if(points[index].index == vertex) {
            break;
        }
    }

    /* check if the segment is found */
    for(i = 0; i < segments.size(); i++) {

        if(points[segments[i][0]].index == vertex) {
            for(j = 0; j < size; j++) {
                if(points[j].index == points[segments[i][1]].index) {
                    break;
                }
            }
            tmp = new int [2];
            tmp[0] = index;
            tmp[1] = j;
            edges.push_back(tmp);
        }
        else if(points[segments[i][1]].index == vertex) {
            for(j = 0; j < size; j++) {
                if(points[j].index == points[segments[i][0]].index) {
                    break;
                }
            }
            tmp = new int [2];
            tmp[0] = index;
            tmp[1] = j;
            edges.push_back(tmp);
        }
    }
    return edges;
}

/* Finds all disjoint edges between two polygons
 * @param A, the nested polygon
 * @param B, the other polygon
 * @returns disjoint, the vector of disjoint edges
 */
vector<int *> disjoint_edges(Polygon A, Polygon B)
{
    vector<int *> disjoint;
    vector<int *> segments;
    int i = 0;
    int count = 0;
    int *tmp;
    /* create temporary segments vector */
    for(i = 0; i < B.shape.size() - 1; i++) {
        tmp = new int [2];
        tmp[0] = B.shape[i];
        tmp[1] = B.shape[i + 1];
        segments.push_back(tmp);
    }
    for(i = 0; i < A.shape.size() - 1; i++) {
        if(segment_match(segments, A.shape[i], A.shape[i + 1]) == -1) {
            tmp = new int [2];
            tmp[0] = i;
            tmp[1] = i + 1;
            disjoint.push_back(tmp);
        }
    }
    return disjoint;
}

/* Finds all shared edges between two polygons
 * @param A, the first polygon
 * @param B, the second polygon
 * @returns shared, the vector of shared edges
 */
vector<int *> shared_edges(Polygon A, Polygon B)
{
    vector<int *> shared;
    vector<int *> segments;
    int i = 0;
    int count = 0;
    int *tmp;
    /* create temporary segments vector */
    for(i = 0; i < B.shape.size() - 1; i++) {
        tmp = new int [2];
        tmp[0] = B.shape[i];
        tmp[1] = B.shape[i + 1];
        segments.push_back(tmp);
    }
    for(i = 0; i < A.shape.size() - 1; i++) {
        if(segment_match(segments, A.shape[i], A.shape[i + 1]) > -1) {
            tmp = new int [2];
            tmp[0] = i;
            tmp[1] = i + 1;
            shared.push_back(tmp);
        }
    }
    return shared;
}

/* finds all shared points between two polygons
 * @param A, the first polygon
 * @param B, the second polygon
 * @return shared, the vector of shared points
 */
vector<int> shared_points(Polygon A, Polygon B)
{
    vector<int> shared;
    for(int i = 0; i < A.shape.size() - 1; i++) {
        if(shape_match(B.shape, A.shape[i]) > -1) {
            shared.push_back(i);
        }
    }
    return shared;
}

/* uses BFS traversal to reach all edges of a seed, called recursively
 * @param segments, all edges which form graph G
 * @param processed, the nodes that have already been visited completely
 * @param seed, the seed to traverse before processing child nodes
 * @param points, the dataset of points
 * @param size, the size of the dataset
 * @return processed, the vector of processed indices
 */
vector<int> breadth_first_index_search(vector<int *> segments, vector<int> *processed, vector<int> seed, Point *points, int size)
{
    vector<int> indices;
    int i;

    for(i = 0; i < seed.size(); i++) {
        /* processed is updated at the end to mark seed[i] as processed */
        printf("TESTING VERTEX = %d ...\n", seed[i]);
        vector<int> index_buff = index_search(segments, processed, seed[i], points, size);
        printf("NEIGHBORS : ");
        for(int index : index_buff) {
            printf("%d ", points[index].index);
            /* insert if not already in indices */
            if(find(indices.begin(), indices.end(), index) == indices.end()) {
                indices.push_back(index);
            }
        }
        printf("\n");
    }

    /* remove the current seed from the indices to generate a new seed */
    for(int vert : seed) {
        for(i = 0; i < indices.size(); i++) {
            if(points[indices[i]].index == vert) {
                break;
            }
        }
        if(i < indices.size()) {
            indices.erase(indices.begin() + i);
        }
    }

    /* base case */
    if(indices.size() <= 0) {
        return *processed;
    }

    return breadth_first_index_search(segments, processed, indices, points, size);
}

/* searches through a vector of segments for all recursive matching end or
 * @param segments, all edges which form graph G
 * @param processed, the nodes that have already been visited completely
 * @param vertex, the point.index value to start from (not an index)
 * @param points, the dataset of points
 * @param size, the size of the dataset
 * @return indices, all the indices found to have a unique match
 */
vector<int> index_search(vector<int *> segments, vector<int> *processed, int vertex, Point *points, int size)
{
    vector<int> indices;

    if(find(processed->begin(), processed->end(), vertex) == processed->end()) {
        /* find all segements with an endpoint at the vertex */
        for(int i = 0; i < segments.size(); i++) {
            if(points[segments[i][0]].index == vertex) {
                /* don't process if the segment connects to a node already been visited */
                if(find(processed->begin(), processed->end(), points[segments[i][1]].index) == processed->end()) {
                    indices.push_back(segments[i][1]);
                }
            }
            else if(points[segments[i][1]].index == vertex) {
                /* don't process if the segment connects to a node already been visited */
                if(find(processed->begin(), processed->end(), points[segments[i][0]].index) == processed->end()) {
                    indices.push_back(segments[i][0]);
                }
            }
        }

        /* mark the vertex as processed */
        processed->push_back(vertex);
    }

    return indices;
}

/* check if two vectors are intersecting
 * @param V1, the first vector
 * @param V2, the second vector
 * @return true if the vectors intersect, false otherwise
 */
bool intersection(Vector V1, Vector V2)
{
    /* return false if any of the points are equal */
    if(V1.start.equals(V2.start) || V1.start.equals(V2.end) || V1.end.equals(V2.start) || V1.end.equals(V2.end)) {
        return false;
    }

    /* sort the points of V1 */
    if(V1.i == 0) { // use y values
        if(V1.end.y < V1.start.y) {
            V1 = Vector(V1.name, V1.end, V1.start);
        }
    }
    else { // use x values
        if(V1.end.x < V1.start.x) {
            V1 = Vector(V1.name, V1.end, V1.start);
        }
    }

    /* sort the points of V2 */
    if(V2.i == 0) { // use y values
        if(V2.end.y < V2.start.y) {
            V2 = Vector(V2.name, V2.end, V2.start);
        }
    }
    else { // use x values
        if(V2.end.x < V2.start.x) {
            V2 = Vector(V2.name, V2.end, V2.start);
        }
    }

    /* parallel test */
    if(determinant(V1, V2) != 0) {
        double y = 0.0;
        double x = 0.0;
        double m1 = 0.0;
        double m2 = 0.0;
        double b1 = 0.0;
        double b2 = 0.0;
        /* find which intercept to use, x or y */
        if(V1.i == 0) { // use x intercept
            /* create first equation */
            y = V1.start.x;
            m1 = (V1.end.x - V1.start.x) / (V1.end.y - V1.start.y);
            x = V1.start.y;
            b1 = y - m1 * x;
        }
        else { // use y intercept
            /* create second equation */
            y = V1.start.y;
            m1 = (V1.end.y - V1.start.y) / (V1.end.x - V1.start.x);
            x = V1.start.x;
            b1 = y - m1 * x;
        }
        if(V2.i == 0) { // use x intercept
            /* create first equation */
            y = V2.start.x;
            m2 = (V2.end.x - V2.start.x) / (V2.end.y - V2.start.y);
            x = V2.start.y;
            b2 = y - m2 * x;
        }
        else { // use y intercept
            /* create second equation */
            y = V2.start.y;
            m2 = (V2.end.y - V2.start.y) / (V2.end.x - V2.start.x);
            x = V2.start.x;
            b2 = y - m2 * x;
        }
        /* solve linear system */
        x = (b2 - b1) / (m1 - m2);
        y = m1 * x + b1;

        /* make final tests */
        Point t = Point(x, y, -1);
        Vector T = Vector("T", V1.start, t);
        if(same_direction(T, V1) && (T.length < V1.length)) { // vector T is inside of vector V1
            T = Vector("T", V2.start, t);
            if(same_direction(T, V2) && (T.length < V2.length)) { // vector T is also inside of vector V2
                return true;
            }
        }
    }

    return false;
}

/* tests whether two vectors are going in the same direction or not
 * @param V1, the first vector
 * @param V2, the second vector
 * @return true if they are in the same direction, false otherwise
 */
bool same_direction(Vector V1, Vector V2) {
    if((V1.i == V2.i) || (V1.i < 0 && V2.i < 0) || (V1.i > 0 && V2.i > 0)) {
        if((V1.j == V2.j) || (V1.j < 0 && V2.j < 0) || (V1.j > 0 && V2.j > 0)) {
            return true;
        }
    }

    return false;
}

/* adds two polygons together to create a single polygon
 * @param A, the first polygon
 * @param B, the second polygon
 * @param points, the point array of datapoints
 * @return C, the created polygon
 */
Polygon add_polygons(Polygon A, Polygon B, Point *points)
{
    Polygon C; //A + B
    int i = 0;
    int j = 0;
    int a1 = 0;
    int a2 = 0;
    int b1 = 0;
    int b2 = 0;
    int tmp = 0;

    /* print current shapes */
    printf("\nBEFORE: A = %0.2lf, B = %0.2lf\n", A.perimeter, B.perimeter);
    C.perimeter = 0.0;
    for(i = 0; i < A.shape.size() - 1; i++) {
        /* find the segment that matches */
        if(((b1 = shape_match(B.shape, A.shape[i])) > -1) && ((b2 = shape_match(B.shape, A.shape[i + 1])) > -1)) {
            printf("b1 = %d, b2 = %d\n", b1, b2);
            a1 = i;
            a2 = i + 1;
            if(a1 == A.shape.size() - 1) {
                a1 = 0;
            }
            else if(a2 == A.shape.size() - 1) {
                a2 = 0;
            }
            printf("a1 = %d, a2 = %d\n", a1, a2);
            /* check which direction to traverse */
            if(b1 == ((b2 + 1) % (B.shape.size() - 1))) {
                /* traverse A to the right */
                if(a1 == (a2 + 1) % (A.shape.size() - 1)) {
                    printf("A: RIGHT\n");
                    i = a1;
                    while(i != a2) {
                        C.shape.push_back(A.shape[i]);
                        C.perimeter += distance_p(points[A.shape[i]], points[A.shape[i + 1]]);
                        printf("%d->", points[A.shape[i]].index);
                        i++;
                        if(i > A.shape.size() - 2)
                            i = 0;
                    }
                }
                /* traverse A to the left */
                else {
                    printf("A: LEFT\n");
                    i = a1;
                    while(i != a2) {
                        C.shape.push_back(A.shape[i]);
                        C.perimeter += distance_p(points[A.shape[i]], points[A.shape[i + 1]]);
                        printf("%d->", points[A.shape[i]].index);
                        i--;
                        if(i < 0)
                            i = A.shape.size() - 2;
                    }
                }
                /* traverse B to the left */
                printf("B: LEFT\n");
                j = b2;
                while(j != b1) {
                    C.shape.push_back(B.shape[j]);
                    C.perimeter += distance_p(points[B.shape[j]], points[B.shape[j + 1]]);
                    printf("%d->", points[B.shape[j]].index);
                    j--;
                    if(j < 0)
                        j = B.shape.size() - 2;
                }
                C.shape.push_back(B.shape[j]);
                C.perimeter += distance_p(points[B.shape[j]], points[B.shape[j + 1]]);
                printf("%d\n", points[B.shape[j]].index);
            }
            else {
                /* traverse A to the right */
                if(a1 == (a2 + 1) % (A.shape.size() - 1)) {
                    printf("A: RIGHT\n");
                    i = a1;
                    while(i != a2) {
                        C.shape.push_back(A.shape[i]);
                        C.perimeter += distance_p(points[A.shape[i]], points[A.shape[i + 1]]);
                        printf("%d->", points[A.shape[i]].index);
                        i++;
                        if(i > A.shape.size() - 2)
                            i = 0;
                    }
                }
                /* traverse A to the left */
                else {
                    printf("A: LEFT\n");
                    i = a1;
                    while(i != a2) {
                        C.shape.push_back(A.shape[i]);
                        C.perimeter += distance_p(points[A.shape[i]], points[A.shape[i + 1]]);
                        printf("%d->", points[A.shape[i]].index);
                        i--;
                        if(i < 0)
                            i = A.shape.size() - 2;
                    }
                }
                /* traverse B to the right */
                printf("B: RIGHT\n");
                j = b2;
                while(j != b1) {
                    C.shape.push_back(B.shape[j]);
                    C.perimeter += distance_p(points[B.shape[j]], points[B.shape[j + 1]]);
                    printf("%d->", points[B.shape[j]].index);
                    j++;
                    if(j > B.shape.size() - 2)
                        j = 0;
                }
                C.shape.push_back(B.shape[j]);
                C.perimeter += distance_p(points[B.shape[j]], points[B.shape[j + 1]]);
                printf("%d\n", points[B.shape[j]].index);
            }
            break;
        }
    }

    /* print altered shapes */
    printf("\nNEW SHAPE: C = %0.2lf\n", C.perimeter);

    return C;
}

/* removes whichever polygon is smaller
 * @param A, the first polygon
 * @param B, the second polygon
 * @param points, the point array of datapoints
 * @return C, the created polygon
 */
Polygon sub_polygons(Polygon A, Polygon B, Point *points)
{
    Polygon C;
    Polygon tmp_polygon;
    vector<int *> found_edges;
    vector<int *> segments;
    int *tmp;
    int i = 0;
    int a1 = 0;
    int a2 = 0;
    int b1 = 0;
    int b2 = 0;
    int count = 0;
    /* ensure A > B */
    if(A.perimeter < B.perimeter) {
        tmp_polygon = A;
        A = B;
        B = tmp_polygon;
    }
    found_edges = disjoint_edges(B, A);
    printf("\nFound edges: (%zu)\n", found_edges.size());
    for(i = 0; i < found_edges.size(); i++) {
        printf("(%d,%d) = ", found_edges[i][0], found_edges[i][1]);
        printf("<%d,%d>\n", points[B.shape[found_edges[i][0]]].index, points[B.shape[found_edges[i][1]]].index);
    }
    printf("\nPolygon A:\n");
    for(i = 0; i < A.shape.size() - 1; i++) {
        printf("%d->", points[A.shape[i]].index);
    }
    printf("%d\n", points[A.shape[i]].index);
    printf("\nPolygon B:\n");
    for(i = 0; i < B.shape.size() - 1; i++) {
        printf("%d->", points[B.shape[i]].index);
    }
    printf("%d\n\n", points[B.shape[i]].index);
    for(i = 0; i < found_edges.size(); i++) {
        a1 = shape_match(A.shape, B.shape[found_edges[i][0]]);
        a2 = shape_match(A.shape, B.shape[found_edges[i][1]]);
        b1 = found_edges[i][0];
        b2 = found_edges[i][1];
        printf("b1 = %d, b2 = %d\n", b1, b2);
        printf("a1 = %d, a2 = %d\n", a1, a2);
        /* find which way to traverse */
        if(b1 == ((b2 + 1) % (B.shape.size() - 1))) {
            if(a1 < a2) {
                /* traverse left */
                printf("SUB A LEFT: ");
                i = a1;
                while(i != a2) {
                    C.shape.push_back(A.shape[i]);
                    C.perimeter += distance_p(points[A.shape[i]], points[A.shape[i + 1]]);
                    printf("%d->", points[A.shape[i]].index);
                    i--;
                    if(i < 0) {
                        i = A.shape.size() - 2;
                    }
                }
            }
            else {
                /* traverse right */
                printf("SUB A RIGHT: ");
                i = a1;
                while(i != a2) {
                    C.shape.push_back(A.shape[i]);
                    C.perimeter += distance_p(points[A.shape[i]], points[A.shape[i + 1]]);
                    printf("%d->", points[A.shape[i]].index);
                    i++;
                    if(i > A.shape.size() - 2) {
                        i = 0;
                    }
                }
            }
            /* traverse right */
            printf("SUB B RIGHT: ");
            i = b1;
            while(i != b2) {
                C.shape.push_back(B.shape[i]);
                C.perimeter += distance_p(points[B.shape[i]], points[B.shape[i + 1]]);
                printf("%d->", points[B.shape[i]].index);
                i++;
                if(i > B.shape.size() - 2) {
                    i = 0;
                }
            }
            C.shape.push_back(B.shape[b2]);
            C.perimeter += distance_p(points[B.shape[b2]], points[B.shape[b1]]);
            printf("%d\n", points[B.shape[i]].index);
        }
        else {
            if(a1 < a2) {
                /* traverse left */
                printf("SUB A LEFT: ");
                i = a1;
                while(i != a2) {
                    C.shape.push_back(A.shape[i]);
                    C.perimeter += distance_p(points[A.shape[i]], points[A.shape[i + 1]]);
                    printf("%d->", points[A.shape[i]].index);
                    i--;
                    if(i < 0) {
                        i = A.shape.size() - 2;
                    }
                }
            }
            else {
                /* traverse right */
                printf("SUB A RIGHT: ");
                i = a1;
                while(i != a2) {
                    C.shape.push_back(A.shape[i]);
                    C.perimeter += distance_p(points[A.shape[i]], points[A.shape[i + 1]]);
                    printf("%d->", points[A.shape[i]].index);
                    i++;
                    if(i > A.shape.size() - 2) {
                        i = 0;
                    }
                }
            }
            printf("SUB B LEFT: ");
            i = b2;
            while(i != b1) {
                C.shape.push_back(B.shape[i]);
                C.perimeter += distance_p(points[B.shape[i]], points[B.shape[i + 1]]);
                printf("%d->", points[B.shape[i]].index);
                i--;
                if(i < 0) {
                    i = B.shape.size() - 2;
                }
            }
            C.shape.push_back(B.shape[b1]);
            C.perimeter += distance_p(points[B.shape[b1]], points[B.shape[b2]]);
        }
        printf("\n");
    }
    return C;
}

/* finds the angle between another vector
 * @param V1, the first vector
 * @param V2, the second vector
 * @return the angle in radians
 */
double angle(Vector V1, Vector V2) {
    return (acos(dot_product(V1, V2) / (V1.length * V2.length)));
}

/* calculates distance given two points
 * @param P1, the start point
 * @param P2, the end point
 * @return the distance between the points
 */
double distance_p(Point P1, Point P2) {
    return sqrt(pow(P2.x - P1.x, 2) + pow(P2.y - P1.y, 2));
}

/* finds the distance between the endpoints of two vectors
 * @param V1, the first vector
 * @param V2, the second vector
 * @return the distance between the endpoints of the vectors
 */
double distance_v(Vector V1, Vector V2) {
    return sqrt(pow(V2.i - V1.i, 2) + pow(V2.j - V1.j, 2));
}

/* finds the dot product of two vectors
 * @param V1, the first vector
 * @param V2, the second vector
 * @return the dot product of the vectors
 */
double dot_product(Vector V1, Vector V2) {
    return ((V1.i * V2.i) + (V1.j * V2.j));
}

/* finds the projection_V1(V2) of two vectors
 * @param V1, the first vector
 * @param V2, the second vector
 * @return projection_V1(V2)
 */
Vector projection(Vector V1, Vector V2) {
    double coeff = (dot_product(V1, V2) / (V1.length * V1.length));
    double i = V1.i * coeff;
    double j = V1.j * coeff;
    Vector result = Vector(V1);
    result.end = Point(result.start);
    result.end.offset(i, j);
    result.refresh();
    return result;
}

/* calculates the ciruclar gradient for a point p
 * @param center, the center point of the circle
 * @param radius, the radius of the circle
 * @param p, the point to use in calcualtion
 * @return grad, the gradient vector
 */
Vector circular_gradient(Point center, double radius, Point p) {
    /* Using the equation f(x,y) = [(x - x_0)^2 + (y - y_0)^2]/r^2
     * grad(f) = d/dx[(x - x_0)^2 + (y - y_0)^2)]i/r^2 +
     *           d/dy[(x - x_0)^2 + (y - y_0)^2]j/r^2
     *         = d/dx[x^2 - 2xx_0 + x_0^2 + y^2 - 2yy_0 + y_0^2]i/r^2 +
     *           d/dy[x^2 - 2xx_0 + x_0^2 + y^2 - 2yy_0 + y_0^2]j/r^2
     *         = (2x - 2x_0)i/r^2 + (2y - 2y_0)j/r^2 */
    Vector grad = Vector("G", center, center);
    grad.i = (2 * (p.x) - 2 * center.x) / pow(radius, 2);
    grad.j = (2 * (p.y) - 2 * center.y) / pow(radius, 2);
    grad.end.offset(grad.i, grad.j);
    grad.refresh();
    return grad;
}

/* finds the intersection point between two vectors
 * @param V1, one of the vectors to test
 * @param V2, the other vector to test
 * @return m, a point with index 0 if found, -1 if not found
 */
Point find_intersection(Vector V1, Vector V2) {
    double y = 0.0;
    double x = 0.0;
    double m1 = 0.0;
    double m2 = 0.0;
    double b1 = 0.0;
    double b2 = 0.0;

    /* find intercept point */
    if(V1.i == 0) { // shift axis to V1
        /* create equation for V2 */
        if(V2.i == 0) {
            return Point(-1, -1, -1); // invalid state
        }
        else {
            y = V2.start.y;
            m2 = (V2.end.y - V2.start.y) / (V2.end.x - V2.start.x);
            x = V2.start.x;
            b2 = y - m2 * x;
        }

        x = V1.start.x;
        y = m2 * x + b2;
    }
    else if(V2.i == 0) { // shift axis to V2
        /* create equation for V2 */
        y = V1.start.y;
        m2 = (V1.end.y - V1.start.y) / (V1.end.x - V1.start.x);
        x = V1.start.x;
        b2 = y - m2 * x;

        x = V2.start.x;
        y = m2 * x + b2;
    }
    else { // solve regularly
        /* create equation for V1 */
        y = V1.start.y;
        m1 = (V1.end.y - V1.start.y) / (V1.end.x - V1.start.x);
        x = V1.start.x;
        b1 = y - m1 * x;

        /* create equation for V2 */
        y = V2.start.y;
        m2 = (V2.end.y - V2.start.y) / (V2.end.x - V2.start.x);
        x = V2.start.x;
        b2 = y - m2 * x;

        x = (b2 - b1) / (m1 - m2);
        y = m1 * x + b1;
    }

    /* check if vectors intersects */
    Point m = Point(x, y, 0);
    Vector M = Vector("M", V1.start, m);
    if(same_direction(M, V1) && (M.length < V1.length)) { // vector M is inside of vector V1
        M = Vector("M", V2.start, m);
        if(same_direction(M, V2) && (M.length < V2.length)) { // vector M is also inside of vector V2
            /* the vectors intersect, return point */
            return m;
        }
    }

    return Point(-1, -1, -1);
}

/* returns the determinant between two vectors
 * @param V1, the first vector
 * @param V2, the second vector
 * @return the determinant, i_1*j_2 - j_1*i_2
 */
double determinant(Vector V1, Vector V2) {
    return ((V1.i * V2.j) - (V2.i * V1.j));
}

/* detects if two segments overlap
 * @param V1, the first segment to test
 * @param V2, the second segment to test
 * @return true if they overlap, false otherwise
 */
bool overlap(Vector V1, Vector V2) {
    /* sort the points of V1 */
    if(V1.i == 0) { // use y values
        if(V1.end.y < V1.start.y) {
            V1 = Vector(V1.name, V1.end, V1.start);
        }
    }
    else { // use x values
        if(V1.end.x < V1.start.x) {
            V1 = Vector(V1.name, V1.end, V1.start);
        }
    }

    /* sort the points of V2 */
    if(V2.i == 0) { // use y values
        if(V2.end.y < V2.start.y) {
            V2 = Vector(V2.name, V2.end, V2.start);
        }
    }
    else { // use x values
        if(V2.end.x < V2.start.x) {
            V2 = Vector(V2.name, V2.end, V2.start);
        }
    }

    /* parallel test */
    if(determinant(V1, V2) == 0) {
        double y = 0.0;
        double x = 0.0;
        double m1 = 0.0;
        double m2 = 0.0;
        double b1 = 0.0;
        double b2 = 0.0;
        /* find which intercept to use, x or y */
        if(V1.i == 0) { // use x intercept
            /* create first equation */
            y = V1.start.x;
            m1 = (V1.end.x - V1.start.x) / (V1.end.y - V1.start.y);
            x = V1.start.y;
            b1 = y - m1 * x;
            /* create second equation */
            y = V2.start.x;
            m2 = (V2.end.x - V2.start.x) / (V2.end.y - V2.start.y);
            x = V2.start.y;
            b2 = y - m2 * x;
        }
        else { // use y intercept
            /* create first equation */
            y = V1.start.y;
            m1 = (V1.end.y - V1.start.y) / (V1.end.x - V1.start.x);
            x = V1.start.x;
            b1 = y - m1 * x;
            /* create second equation */
            y = V2.start.y;
            m2 = (V2.end.y - V2.start.y) / (V2.end.x - V2.start.x);
            x = V2.start.x;
            b2 = y - m2 * x;
        }
        if(b1 == b2) {
            if(V1.i == 0) { // use y values
                if(V1.start.y <= V2.start.y) { // P1 <= P3
                    if(V2.start.y < V1.end.y) { // P1 --> P3 --> P2 --> P4
                        return true;
                    }
                }
                else if(V2.start.y <= V1.start.y) { // P3 <= P1
                    if(V1.start.y < V2.end.y) { // P3 --> P1 --> P4 --> P2
                        return true;
                    }
                }
            }
            else { // use x values
                if(V1.start.x <= V2.start.x) { // P1 < P3
                    if(V2.start.x < V1.end.x) { // P1 --> P3 --> P2 --> P4
                        return true;
                    }
                }
                else if(V2.start.x <= V1.start.x) { // P3 <= P1
                    if(V1.start.x < V2.end.x) { // P3 --> P1 --> P4 --> P2
                        return true;
                    }
                }
            }
        }
    }

    return false;
}

/* finds the convex hull of a given set of points using Graham's scan algorithm
 * @param points, the set of datapoints
 * @param size, the size of the points array
 * @return polygon, the polygon representing the convex hull
 */
Polygon find_convex_hull(Point *points, int size) {
    int i;
    int j;

    /* find the lowest, right-most point */
    double min_y = DBL_MAX;
    double max_x = DBL_MIN;
    int k = -1;
    for(i = 0; i < size; i++) {
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
    for(i = 0; i < size; i++) {
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
    shape.push_back(k);
    shape.push_back(point_match(points, size, M.end.index));

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
        shape.push_back(point_match(points, size, M.end.index));
    }

    /* return polygon */
    Polygon polygon = Polygon(shape, points);
    return polygon;
}

/* returns all wall segments, including crosses
 * @param points, the point array of datapoints
 * @param size, the size of the array
 * @return the vector of segments calculated
 */
vector<int *> all_w_segments(Point *points, int size) {

    /////////////////////////////////
    // INITIALIZATION OF VARIABLES //
    /////////////////////////////////

    Polygon convex_hull = find_convex_hull(points, size);
    vector<int *> w_segments = convex_hull.segments; // output segments
    vector<int *> t_segments; // test segments
    int *tmp_segment = new int [2];
    int i = 0;
    int j = 0;
    int k = 0;

    //////////////////////////
    // CALCULATE W-SEGMENTS //
    //////////////////////////

    /* create the interval based off the smallest segment
    double interval = DBL_MAX;
    for(i = 0; i < size; i++) {
        for(j = 0; j < size; j++) {
            if(i == j) {
                continue;
            }
            if(distance_p(points[i], points[j]) < interval) {
                interval = distance_p(points[i], points[j]);
            }
        }
    }

    interval /= 5;*/

    /* test all line segments that are inside the base but not a part of the base */
    for(i = 0; i < size; i++) {
        for(j = 0; j < size; j++) {
            if(i == j) {
                continue;
            }
            if(segment_match(w_segments, i, j) == -1) {
                /* record segment */
                tmp_segment = new int [2];
                tmp_segment[0] = i;
                tmp_segment[1] = j;
                t_segments.push_back(tmp_segment);
            }
        }
    }

    /* bubble sort test segments by smallest length */
   for(i = 0; i < t_segments.size(); i++) {
        for(j = t_segments.size() - 1; j > i; j--) {
            if(distance_p(points[t_segments[j][0]], points[t_segments[j][1]]) < distance_p(points[t_segments[j - 1][0]], points[t_segments[j - 1][1]])) {
                tmp_segment = new int [2];
                tmp_segment[0] = t_segments[j][0];
                tmp_segment[1] = t_segments[j][1];
                t_segments[j] = t_segments[j - 1];
                t_segments[j - 1] = tmp_segment;
            }
        }
    } 

    for(i = 0; i < t_segments.size(); i++) {
        /*if((t_segments[i][0] == 3 && t_segments[i][1] == 5) || (t_segments[i][0] == 5 && t_segments[i][1] == 3))
            break;*/
        Vector L = Vector("L", points[t_segments[i][0]], points[t_segments[i][1]]);

        /* normalize L vector */
        if(L.i == 0) { // use y values
            if(L.start.y > L.end.y) {
                L = Vector(L.name, L.end, L.start);
            }
        }
        else { // use x values
            if(L.start.x > L.end.x) {
                L = Vector(L.name, L.end, L.start);
            }
        }

        if(segment_match(w_segments, t_segments[i][0], t_segments[i][1]) == -1) {
            /* if the line is valid check for crosses */
            if(test_w_segment(w_segments, L, 0, points, size)) { // interval value is not necessary...
                /* record segment */
                tmp_segment = new int [2];
                tmp_segment[0] = t_segments[i][0];
                tmp_segment[1] = t_segments[i][1];

                w_segments.push_back(tmp_segment);
            }
        }
    }

    return w_segments;
}

/* creates a normalized bijection vector to a
 * @param L, the tested segment
 * @param points, the set of datapoints
 * @param n, the size of the points array
 * @return true if the segment is validated, false otherwise
 */
bool test_w_segment(vector<int *> segments, Vector L, double interval, Point *points, int n) {
    Vector V1;
    Vector V2;
    vector<int> top; // stores the left-side indices
    vector<int> bottom; // stores the right-side indices
    int *ignore = new int [n]; // ignore array
    int i = 0;

    for(i = 0; i < n; i++) {
        ignore[i] = -1;
    }

    /* find the midpoint of L */
    Point m = Point(L.start.x, L.start.y, -1);
    m.offset(L.i / 2, L.j / 2);

    /* store radius */
    double radius = distance_p(L.start, m);

    /* run tests on point set */
    for(i = 0; i < n; i++) {
        Point p = Point(points[i]);

        /* skip starting at L.start and L.end */
        if(p.equals(L.start) || p.equals(L.end)) {
            continue;
        }

        /* store point if it is within a hemisphere */
        if(distance_p(p, m) <= radius) {
            Vector P = Vector("P", L.start, p);

            if(determinant(L, P) > 0) {
                top.push_back(i);
            }
            else if(determinant(L, P) < 0) {
                bottom.push_back(i);
            }
            else {
                return false; // invalid, point is on L
            }

            /* iterate through test segments */
            if(top.size() > 0 && bottom.size() > 0) {
                for(int t : top) {
                    for(int b : bottom) {
                        /* check if straddling point pairs have not been checked */
                        if(ignore[t] != b) {
                            Vector Q = Vector("Q", points[t], points[b]);
                            vector<int> *processed;

                            /* calculate degree of Q */
                            processed = new vector<int>;
                            vector<int> edges_q1 = index_search(segments, processed, Q.start.index, points, n);
                            vector<int> edges_q2 = index_search(segments, processed, Q.end.index, points, n);
                            int degree_q = edges_q1.size() + edges_q2.size();

                            /* calculate degree of L */
                            processed = new vector<int>;
                            vector<int> edges_l1 = index_search(segments, processed, L.start.index, points, n);
                            vector<int> edges_l2 = index_search(segments, processed, L.end.index, points, n);
                            int degree_l = edges_l1.size() + edges_l2.size();

                            /* L is invalid if degree of Q > degree of L */
                            if(degree_q > degree_l) {
                                printf("INVALID!\n");
                                return false; // invalid
                            }

                            /* add the points as ignored */
                            ignore[t] = b;
                            ignore[b] = t;
                        }
                    }
                }
            }
        }
        else {
            //printf("NOT IN BIJECTION RANGE!\n");
        }

    }

    return true; //valid
}

/* generates all shortest paths from W
 * @param segments, all edges which form graph G
 * @param points, the dataset of points
 * @param n, the size of the dataset
 * @return paths, the vector of shortest paths
 */
vector<Polygon> generate_final_paths(vector<int *> segments, Point *points, int n) {
    vector<Polygon> polygons; // all polygons created
    Polygon path; // current path
    unordered_map<string, int> map; // map of processed segment to polygon index
    unordered_map<string, vector<string>> cross; // map of segment to crossing value
    string key; // key used for storing maps
    string value; // value used for storing maps
    ostringstream stream; // formatting for key
    vector<int *> seeds; // set of segments to test
    vector<int *> set; // working set of segments during test
    int i = 0;
    int j = 0;
    int k = 0;

    /* store matrix of crossing segments */
    if(segments.size() % 2 == 0) {
        k = segments.size() / 2;
    }
    else {
        k = (segments.size() / 2) + 1;
    }
    for(int *segment_a : segments) {
        Vector A = Vector("A", points[segment_a[0]], points[segment_a[1]]);

        for(int *segment_b : segments) {
            /* skip if any of the points are equal */
            if(segment_b[0] == segment_a[0] || segment_b[1] == segment_a[1] || segment_b[0] == segment_a[1] || segment_b[1] == segment_a[0]) {
                continue;
            }
            Vector B = Vector("B", points[segment_b[0]], points[segment_b[1]]);
            Point t = find_intersection(A, B);
            if(t.index == 0) {
                /* create key */
                stream.flush();
                if(segment_a[0] < segment_a[1]) {
                    stream << segment_a[0] << "," << segment_a[1];
                }
                else {
                    stream << segment_a[1] << "," << segment_a[0];
                }
                key = stream.str();

                /* create value */
                stream.flush();
                if(segment_b[0] < segment_b[1]) {
                    stream << segment_b[0] << "," << segment_b[1];
                }
                else {
                    stream << segment_b[1] << "," << segment_b[0];
                }
                value = stream.str();

                /* store data */
                vector<string> values;
                unordered_map<string,vector<string>>::const_iterator item;

                item = cross.find(key);
                if (item != cross.end()) {
                    values = cross.at(key);
                }
                values.push_back(value);
                cross.insert({key, values});

                item = cross.find(value);
                if (item != cross.end()) {
                    values = cross.at(value);
                }
                values.push_back(key);
                cross.insert({value, values});
            }
        }
    }

    for(int *segment : segments) {
        vector<string> values;
        unordered_map<string,vector<string>>::const_iterator item;

        /* create key */
        stream.flush();
        if(segment[0] < segment[1]) {
            stream << segment[0] << "," << segment[1];
        }
        else {
            stream << segment[1] << "," << segment[0];
        }
        key = stream.str();

        /* check if there is a cross */
        item = cross.find(key);
        if(item != cross.end()) {
            /* found cross */
        }
    }

    set = segments;

    for(int *segment : set) {
        Polygon smallest = dijkstra_polygon(segments, segment, points, n);
        if(segment[0] < segment[1]) {
            stream << segment[0] << "," << segment[1];
        }
        else {
            stream << segment[1] << "," << segment[0];
        }
        key = stream.str();

        /* check how the polygon should be stored */
        for(i = 0; i < polygons.size(); i++) {
            if(polygons[i].equals(smallest)) {
                map.insert({key, i});
                break;
            }
        }
        if(i == polygons.size()) {
            map.insert({key, polygons.size()});
            polygons.push_back(smallest);
        }

        /* update path */

    }

    return polygons;
}

/* implements a special version of Dijkstra's Algorithm
 * @param segments, all edges which form graph G
 * @param segment, the segment to use as a seed
 * @param points, the dataset of points
 * @param n, the size of the dataset
 */
Polygon dijkstra_polygon(vector<int *> segments, int *segment, Point *points, int n) {
    double *dist = new double [n]; // smallest distances
    vector<int> path; // stores the path of the polygon
    vector<int> *prev = new vector<int> [n]; // vector of sub paths
    vector<int> set; // working set of points
    double smallest = DBL_MAX;
    int index = 0;
    int i = 0;
    int j = 0;

    for(i = 0; i < n; i++) {
        set.push_back(i);
        dist[i] = DBL_MAX;
    }
    dist[segment[0]] = 0; // source to source is 0

    while(set.size() > 0) {
        /* choose the point with the smallest dist value */
        for(j = 0; j < n; j++) {
            if(dist[j] < smallest && find(set.begin(), set.end(), j) != set.end()) {
                smallest = dist[j];
                index = j;
            }
        }

        /* remove points[index] from set of points */
        for(j = 0; j < n; j++) {
            if(set[j] == index) {
                set.erase(set.begin() + j);
                break;
            }
        }

        /* iterate over all neighbors of p */
        vector<int *> edges = edge_search(segments, points[index].index, points, n);
        for(int *edge : edges) {
            /* skip if the prev is going directly from start to end */
            if(edge[0] == segment[0] && edge[1] == segment[1]) {
                continue;
            }

            double alt = dist[index] + distance_p(points[index], points[edge[1]]);
            if(alt < dist[edge[1]]) {
                dist[edge[1]] = alt;
                prev[edge[1]].push_back(index);
            }
        }
        smallest = DBL_MAX;
        index = 0;
    }

    /* store the shortest prev from start to end */
    i = segment[1];
    path.push_back(i);
    while(i != segment[0]) {
        for(j = 0; j < prev[i].size(); j++) {
            path.push_back(prev[i][j]); // store the index in the path
            //printf("PUSHED %d\n", points[prev[i][j]].index);
        }
        i = prev[i][j - 1];
    }

    Polygon polygon = Polygon(path, points);

    return polygon;
}
