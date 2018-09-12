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

/* searches through a vector of segments for all matching end or
 * beginning segments, returning the first and last indices at which
 * it was found in a vector of 2D arrays */
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

/* searches through a vector of polygons for a matching vertex */
int polygons_search(vector<Polygon> polygons, int vertex)
{
    int d = 0;
    /* check if the vertex is found */
    for(int i = 0; i < polygons.size(); i++) {
        if(find(polygons[i].shape.begin(), polygons[i].shape.end(), vertex) != polygons[i].shape.end()) {
            return i;
        }
    }
    return -1;
}

/* searches through a vector of segments for a matching segment */
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

/* searches through a vector of segments for a matching segment */
int index_match(vector<int *> segments, int beginning, int end) {
    for(int i = 0; i < segments.size(); i++) {
        if((segments[i][0] == beginning && segments[i][1] == end) || (segments[i][1] == beginning && segments[i][0] == end)) {
            return i;
        }
    }
    return -1;
}

/* searches through a shape for a duplicate */
int duplicate_search(vector<int> shape)
{
    for(int i = 1; i < shape.size(); i++) {
        /* check if the vertex has a duplicate */
        if(find(shape.begin() + i + 1, shape.end(), shape[i]) != shape.end()) {
            return 1;
        }
    }
    return 0;
}

/* returns true if the vectors intersect, false otherwise */
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

Polygon find_shortest_path(vector<Polygon> polygons, Point *points, int size)
{
    if(polygons.size() == 0) {
        Polygon p;
        return p;
    }
    vector<Polygon> original = polygons;
    vector<Polygon> base;
    vector<Polygon> bridge;
    Polygon tmp_polygon;
    Polygon shortest_path;
    double epsilon = 0.000001;
    double min = DBL_MAX;
    vector<int *> found_edges;
    vector<int *> ref_edges;
    vector<int *> bridge_edges;
    vector<int> found_points;
    int *visited = new int [size] ();
    int *tmp;
    int i = 0;
    int j = 0;
    int k = 0;
    int n = 0;
    int m = 0;
    int accepted = 0;
    int complete = 0;
    int found = 0;
    int sum_edges = 0;
    int sum_points = 0;

    /* bubble sort polygons by perimeter */
    for(i = 0; i < polygons.size(); i++) {
        for(j = polygons.size() - 1; j > i; j--) {
            if(polygons[j].perimeter > polygons[j - 1].perimeter) {
                tmp_polygon = polygons[j];
                polygons[j] = polygons[j - 1];
                polygons[j - 1] = tmp_polygon;
            }
        }
    }
    base.push_back(polygons[0]);
    visit_polygon(visited, polygons[0], points);
    printf("\nADD BASE: ");
    for(k = 0; k < polygons[0].shape.size() - 1; k++) {
        printf("%d->", points[polygons[0].shape[k]].index);
    }
    printf("%d\n", points[polygons[0].shape[k]].index);
    i = 1;
    /* first grow shapes from largest to smallest, without touching the same vertices */
    for(i = 1; i < polygons.size(); i++) {
        /* print bases */
        printf("BASES:\n");
        for(k = 0; k < base.size(); k++) {
            printf("%d: ", k);
            for(n = 0; n < base[k].shape.size() - 1; n++) {
                printf("%d->", points[base[k].shape[n]].index);
            }
            printf("%d\n", points[base[k].shape[n]].index);
        }
        sum_edges = 0;
        sum_points = 0;
        for(j = 0; j < base.size(); j++) {
            found_edges = shared_edges(base[j], polygons[i]);
            sum_edges += found_edges.size();
            found_points = shared_points(base[j], polygons[i]);
            sum_points += found_points.size();
        }
        /* push polygon if there aren't any shared edges or points */
        if((sum_edges == 0) && (sum_points == 0)) {
            printf("\nADD BASE: ");
            for(k = 0; k < polygons[i].shape.size() - 1; k++) {
                printf("%d->", points[polygons[i].shape[k]].index);
            }
            printf("%d\n", points[polygons[i].shape[k]].index);
            base.push_back(polygons[i]);
        }
        /* add polygon if there is only one shared edge */
        else if(sum_edges == 1) {
            /* find which base to add to */
            for(j = 0; j < base.size(); j++) {
                found_edges = shared_edges(base[j], polygons[i]);
                if(found_edges.size() == 1) {
                    break;
                }
            }
            found = 0;
            /*check if any of the new points are contained in another base */
            for(k = 0; k < base.size(); k++) {
                if(k == j) {
                    continue;
                }
                /* loop through all points in each base */
                for(n = 0; n < base[k].shape.size(); n++) {
                    /* loop through all points in the addition */
                    for(m = 0; m < polygons[i].shape.size(); m++) {
                        /* break if a point matches */
                        if((points[base[k].shape[n]].index == points[polygons[j].shape[m]].index) || (points[base[k].shape[n]].index == points[polygons[i].shape[m]].index)) {
                            found = 1;
                            break;
                        }
                    }
                    if(found) {
                        break;
                    }
                }
                if(found) {
                    break;
                }
            }
            /* add shapes if no other base had the new points*/
            if(!found) {
                printf("\nADDITION: ");
                for(k = 0; k < polygons[i].shape.size() - 1; k++) {
                    printf("%d->", points[polygons[i].shape[k]].index);
                }
                printf("%d\n", points[polygons[i].shape[k]].index);
                printf("\nBASE: ");
                for(k = 0; k < base[j].shape.size() - 1; k++) {
                    printf("%d->", points[base[j].shape[k]].index);
                }
                printf("%d\n", points[base[j].shape[k]].index);
                /* add polygon to base */
                base[j] = add_polygons(polygons[i], base[j], points);
                printf("\nNEW BASE: ");
                for(k = 0; k < base[j].shape.size() - 1; k++) {
                    printf("%d->", points[base[j].shape[k]].index);
                }
                printf("%d\n", points[base[j].shape[k]].index);
            }
        }
        /* skip the polygon if there is more than one shared edge */
        else {
            ;
        }
    }
    shortest_path = base[0];
    return shortest_path;
}

/* checks which polygon the segments are contained in */
int accept_polygon(Polygon polygon, vector<int *> segments, Point *points) {
    int i = 0;
    int count = 0;
    printf("Found edges...");
    for(i = 0; i < segments.size(); i++) {
        printf("<%d,%d>", points[segments[i][0]].index, points[segments[i][1]].index);
    }
    printf("\n");
    printf("Checking polygon...");
    for(i = 0; i < polygon.shape.size() - 1; i++) {
        printf("%d->", points[polygon.shape[i]].index);
        if(segment_match(segments, polygon.shape[i], polygon.shape[i + 1]) > -1) {
            count++;
        }
    }
    printf("%d->", points[polygon.shape[i]].index);
    printf("... count = %d, size = %zu\n", count, segments.size());
    if((count > 0) && (count <= segments.size())) {
        return 1;
    }
    return 0;
}

/* finds the smallest neighbour of a polygon */
int smallest_neighbour(vector<Polygon> polygons, Polygon source, int n)
{
    vector<int *> found_edges;
    int i = 0;
    int k = 0;
    int index = 0;
    double min = DBL_MAX;
    for(i = 0; i < polygons.size(); i++) {
        if(i == n) {
            continue;
        }
        found_edges = shared_edges(polygons[i], source);
        if(found_edges.size() > 0) {
            if(min > polygons[i].perimeter) {
                min = polygons[i].perimeter;
                index = i;
            }
        }
    }
    return index;
}

/* A must be the nested polygon */
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

vector<int> shared_points(Polygon A, Polygon B)
{
    vector<int> shared;
    for(int i = 0; i < A.shape.size() - 1; i++) {
        if(shape_search(B.shape, A.shape[i]) > -1) {
            shared.push_back(i);
        }
    }
    return shared;
}

void visit_polygon(int *visited, Polygon polygon, Point *points)
{
    for(int i = 0; i < polygon.shape.size(); i++) {
        visited[points[polygon.shape[i]].index] = 1;
    }
}

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
        if(((b1 = shape_search(B.shape, A.shape[i])) > -1) && ((b2 = shape_search(B.shape, A.shape[i + 1])) > -1)) {
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

/* removes whichever polygon is smaller */
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
        a1 = shape_search(A.shape, B.shape[found_edges[i][0]]);
        a2 = shape_search(A.shape, B.shape[found_edges[i][1]]);
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

/* searches through a shape for a matching vertex */
int shape_search(vector<int> shape, int vertex)
{
    for(int i = 0; i < shape.size() - 1; i++) {
        /* check if the vertex is found */
        if(shape[i] == vertex) {
            return i;
        }
    }
    return -1;
}

/* searches through a shape for a matching edge */
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

/* calculates curvature given structure k */
double calculate_curvature(Vector T1, Vector T2, double tao)
{
    return (distance_v(T1, T2) / angle_t(tao));
}

/* calculates angle given tao */
double angle_t(double tao)
{
    return (acos(tao) + (M_PI / 180));
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

/* prints to the terminal if there is an error assigning memory */
void memory_error(void)
{
    printf("\n\nError assigning memory. Exiting Program. Good Day.\n\n");
}

/* initializes weslean polygon construction
 * @param points, the point array of datapoints
 * @param size, the size of the array
 * @return the vector of polygons calculated
 */
vector<Polygon> init_w_polygons(Point *points, int size) {

    /////////////////////////////////
    // INITIALIZATION OF VARIABLES //
    /////////////////////////////////

    vector<Polygon> polygons;
    Polygon convex_hull = find_convex_hull(points, size);
    vector<int *> w_segments = convex_hull.segments;
    int *tmp_segment = new int [2];
    int i = 0;
    int j = 0;
    int k = 0;

    /////////////////////////////
    // CHECK SIMPLE SHAPE CASE //
    /////////////////////////////

    if(convex_hull.shape.size() - 1 == size) {
        polygons.push_back(convex_hull);
        return polygons;
    }

    //////////////////////////
    // CALCULATE W-SEGMENTS //
    //////////////////////////

    /* find smallest line segment */
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

    /* create the interval based off the smallest segment */
    interval /= 2;
    vector<int *> crosses;

    /* test all line segments that are inside the base but not a part of the base */
    for(i = 0; i < size; i++) {
        for(j = 0; j < size; j++) {
            if(i == j) {
                continue;
            }
            Vector L = Vector("L", points[i], points[j]);
            if(segment_match(w_segments, i, j) == -1) {
                /* if the line is valid check for crosses */
                if(test_w_segment(L, interval, points, size)) {
                    bool push = true;

                    /* check if the segment intersects any of the previously recorded w_segments */
                    for(k = w_segments.size() - 1; k >= 0; k--) {
                        Vector V = Vector("V", points[w_segments[k][0]], points[w_segments[k][1]]);

                        /* check for a non-overlap intersection */
                        if(!overlap(V, L) && intersection(V, L)) {
                            if(segment_match(crosses, w_segments[k][0], w_segments[k][1]) == -1) {
                                crosses.push_back(w_segments[k]);
                            }
                            if(segment_match(crosses, i, j) == -1) {
                                int *tmp_segment = new int [2];
                                tmp_segment[0] = i;
                                tmp_segment[1] = j;
                                crosses.push_back(tmp_segment);
                            }
                            w_segments.erase(w_segments.begin() + k);
                            push = false;
                        }
                    }

                    /* also check if the segment intersects any of the previously deleted w_segments */
                    if(push) {
                        for(k = 0; k < crosses.size(); k++) {
                            Vector V = Vector("V", points[crosses[k][0]], points[crosses[k][1]]);

                            /* always push overlaps, they are handeled later */
                            if(!overlap(V, L) && intersection(V, L)) {
                                if(segment_match(crosses, i, j) == -1) {
                                    int *tmp_segment = new int [2];
                                    tmp_segment[0] = i;
                                    tmp_segment[1] = j;
                                    crosses.push_back(tmp_segment);
                                }
                                push = false;
                                break;
                            }
                        }
                    }

                    /* push the segment if it didn't intersect any other w_segments */
                    if(push) {
                        /* record segment */
                        tmp_segment = new int [2];
                        tmp_segment[0] = i;
                        tmp_segment[1] = j;

                        if(w_segments.size() == 0) {
                            w_segments.push_back(tmp_segment);
                        }
                        else {
                            w_segments = fix_overlap(tmp_segment, w_segments, points);
                        }
                    }
                }
            }
        }
    }

    //////////////////////////
    // CREATE INITIAL HULLS //
    //////////////////////////

    vector<Polygon> w_polygons;
    vector<Polygon> initial_polygons;
    vector<int> processed_points;
    vector<string> processed_hulls;

    /* find the polygon starting at each edge */
    for(int *segment : w_segments) {
        int index;
        vector<int *> edges;

        /* process each index */
        for(i = 0; i < 2; i++) {
            index = points[segment[i]].index;

            /* check if the index is already processed */
            if(find(processed_points.begin(), processed_points.end(), index) == processed_points.end()) {
                processed_points.push_back(index);
                edges = edge_search(w_segments, index, points, size);

                for(int *edge : edges) {
                    /* store created polygons */
                    vector<Polygon> created = create_polygons(edge, w_segments, points, size);

                    /* add each hull to w_polygons if it isn't already added */
                    for(Polygon hull : created) {
                        /* skip the convex hull */
                        if(hull.id == convex_hull.id) {
                            continue;
                        }

                        /* check if we should add the hull to the w_polygons vector */
                        if(find(processed_hulls.begin(), processed_hulls.end(), hull.id) == processed_hulls.end()) {
                            initial_polygons.push_back(hull);
                            processed_hulls.push_back(hull.id);
                        }
                    }
                }
            }
        }
    }

    /////////////////////////////////
    // PROCESS CROSSING W-SEGMENTS //
    /////////////////////////////////

    /* bubble sort crosses for later processing */
    for(i = 0; i < crosses.size(); i++) {
        for(j = crosses.size() - 1; j > i; j--) {
            if(distance_p(points[crosses[j][0]], points[crosses[j][1]]) < distance_p(points[crosses[j - 1][0]], points[crosses[j - 1][1]])) {
                int *tmp_segment = new int [2];
                tmp_segment[0] = crosses[j][0];
                tmp_segment[1] = crosses[j][1];
                crosses[j] = crosses[j - 1];
                crosses[j - 1] = tmp_segment;
            }
        }
    }

    /* recursively generate other w_polygons off initial hulls */
    vector<string> added_hulls;
    for(Polygon hull : initial_polygons) {
        vector<Polygon> created = construct_w_polygons(hull, points, size, processed_hulls, crosses);
        if(created.size() <= 1) {
            if(find(added_hulls.begin(), added_hulls.end(), hull.id) == added_hulls.end()) {
                w_polygons.push_back(hull);
                added_hulls.push_back(hull.id);
            }
        }
        for(Polygon p : created) {
            if(find(processed_hulls.begin(), processed_hulls.end(), p.id) == processed_hulls.end()) {
                w_polygons.push_back(p);
                processed_hulls.push_back(p.id);
            }
        }
    }

    return w_polygons;
}

/* returns all weslean segments, including crosses
 * @param points, the point array of datapoints
 * @param size, the size of the array
 * @return the vector of segments calculated
 */
vector<int *> all_w_segments(Point *points, int size) {

    /////////////////////////////////
    // INITIALIZATION OF VARIABLES //
    /////////////////////////////////

    Polygon convex_hull = find_convex_hull(points, size);
    vector<int *> w_segments = convex_hull.segments;
    int *tmp_segment = new int [2];
    int i = 0;
    int j = 0;
    int k = 0;

    //////////////////////////
    // CALCULATE W-SEGMENTS //
    //////////////////////////

    /* create the interval based off the smallest segment */
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

    interval /= 5;

    /* test all line segments that are inside the base but not a part of the base */
    for(i = 0; i < size; i++) {
        for(j = 0; j < size; j++) {
            if(i == j) {
                continue;
            }
            Vector L = Vector("L", points[i], points[j]);
            if(segment_match(w_segments, i, j) == -1) {
                /* if the line is valid check for crosses */
                if(test_w_segment(L, interval, points, size)) {
                    /* record segment */
                    tmp_segment = new int [2];
                    tmp_segment[0] = i;
                    tmp_segment[1] = j;

                    if(w_segments.size() == 0) {
                        w_segments.push_back(tmp_segment);
                    }
                    else {
                        w_segments = fix_overlap(tmp_segment, w_segments, points);
                    }
                }
            }
        }
    }

    return w_segments;
}

/* constructs weslean polygons given a hull to work off of
 * @param base, the current hull
 * @param points, the point array of datapoints
 * @param size, the size of the array
 * @param processed_hulls, the vector of hulls to ensure uniqueness
 * @param crosses, a vector of all crosses calculated in init
 * @return the vector of polygons recursively calculated
 */
vector<Polygon> construct_w_polygons(Polygon base, Point *points, int size, vector<string> processed_hulls, vector<int *> crosses) {
    int i = 0;
    int j = 0;
    int k = 0;

    printf("CROSSES:\n");
    for(int *cross : crosses) {
        printf("<%d, %d>\n", points[cross[0]].index, points[cross[1]].index);
    }

    vector<Polygon> w_polygons;
    vector<int *> w_segments = base.segments;

    /* add each w_cross in order, and check for w_polygons */
    cout << "BASE: ";
    for(int s : base.shape) {
        cout << points[s].index << " ";
    }
    cout << endl;

    vector<int *> filtered;
    for(i = 0; i < crosses.size(); i++) {
        vector<int *> base_segments = w_segments;
        /* make sure the segment is contained within the polygon */
        if((segment_match(w_segments, crosses[i][0], crosses[i][1]) == -1) && base.contains(points[crosses[i][0]], points, size) && base.contains(points[crosses[i][1]], points, size)) {
            filtered.push_back(crosses[i]);
        }
    }

    cout << "FILTERED: ";
    for(int *s : filtered) {
        cout << "<" << points[s[0]].index << ", " << points[s[1]].index << "> ";
    }
    cout << endl;

    /* check if there is a cross for this base */
    if(filtered.size() == 0) {
        w_polygons.push_back(base);
        return w_polygons;
    }
    if(filtered.size() == 1) {
        vector<Polygon> created = create_polygons(filtered[0], w_segments, points, size);
        if(created.size() <= 1) { // dangling case
            w_polygons.push_back(base);
            return w_polygons;
        }
        else {
            w_polygons.push_back(created[0]);
            w_polygons.push_back(created[1]);
            return w_polygons;
        }
    }

    /* find the crosses and process them */
    for(i = 0; i < crosses.size(); i++) {
        vector<int *> base_segments = w_segments;
        /* make sure the segment is contained within the polygon */
        Point start = points[point_match(points, size, crosses[i][0])];
        if((segment_match(w_segments, crosses[i][0], crosses[i][1]) == -1) && base.contains(points[crosses[i][0]], points, size) && base.contains(points[crosses[i][1]], points, size)) {
            Vector L = Vector("L", points[crosses[i][0]], points[crosses[i][1]]);

            /* fix crosses when encountered */
            printf("CHECKING: <%d, %d>\n", points[crosses[i][0]].index, points[crosses[i][1]].index);
            for(j = 0; j < w_segments.size(); j++) {
                Vector V = Vector("V", points[w_segments[j][0]], points[w_segments[j][1]]);

                /* check for a non-overlap intersection */
                if(!overlap(V, L) && intersection(V, L)) {
                    printf("FIXING CROSS: <%d, %d>, <%d, %d>\n", points[crosses[i][0]].index, points[crosses[i][1]].index, points[w_segments[j][0]].index, points[w_segments[j][1]].index);
                    int *segment = new int [2];
                    int *segment_copy = new int [2];
                    segment[0] = segment_copy[0] = w_segments[j][0];
                    segment[1] = segment_copy[1] = w_segments[j][1];
                    int *cross = new int[2];
                    int *cross_copy = new int[2];
                    cross[0] = cross_copy[0] = crosses[i][0];
                    cross[1] = cross_copy[1] = crosses[i][1];

                    w_segments = base.segments; // reset w_segments
                    w_segments.push_back(cross);
                    vector<Polygon> first;
                    vector<Polygon> second;
                    vector<Polygon> created;
                    Polygon S1;
                    Polygon S2;

                    //TODO: handle dangling crosses

                    /* process first crossing edge for the first polygon */
                    first = create_polygons(cross, w_segments, points, size); // use first segment
                    if(first.size() <= 1) { // dangling cross
                        w_segments = base_segments;
                        break;
                    }

                    /* store the polygon with the greatest perimeter */
                    if(first[0].perimeter > first[1].perimeter) {
                        S1 = first[0];
                    }
                    else {
                        S1 = first[1];
                    }

                    w_segments = base.segments; // reset w_segments
                    w_segments.push_back(segment);

                    /* process second crossing edge for the second polygon */
                    second = create_polygons(segment, w_segments, points, size); // use first segment
                    if(second.size() <= 1) { // dangling cross
                        w_segments = base_segments;
                        break;
                    }

                    /* store the polygon with the greatest perimeter */
                    if(second[0].perimeter > second[1].perimeter) {
                        S2 = second[0];
                    }
                    else {
                        S2 = second[1];
                    }

                    /* keep the polygon with the smallest perimeter */
                    //TODO: DEBUG!
                    if(S1.perimeter < S2.perimeter) {
                        created = first;
                    }
                    else {
                        created = second;
                    }

                    vector<int *> new_crosses;
                    for(k = 0; k < crosses.size(); k++) {
                        if(k != segment_match(crosses, cross_copy[0], cross_copy[1]) && k != segment_match(crosses, segment_copy[0], segment_copy[1])) {
                            new_crosses.push_back(crosses[k]);
                        }
                    }

                    /* add each hull to w_polygons if it isn't already added */
                    for(Polygon hull : created) {
                        /* skip the convex hull */
                        if(hull.id == base.id) {
                            continue;
                        }

                        /* check if we should add the hull to the w_polygons vector */
                        if(find(processed_hulls.begin(), processed_hulls.end(), hull.id) == processed_hulls.end()) {
                            /* recursively generate inner w_polygons */
                            cout << "CREATED HULL: " << hull.id << endl;
                            vector<Polygon> constructed = construct_w_polygons(hull, points, size, processed_hulls, new_crosses);
                            processed_hulls.push_back(hull.id);

                            /* add polygons */
                            vector<Polygon> added;
                            vector<string> added_ids;
                            if(constructed.size() == 0) {
                                added.push_back(hull);
                            }
                            else {
                                for(Polygon p : constructed) {
                                    cout << "NEW... " << p.id << endl;
                                    if(find(added_ids.begin(), added_ids.end(), p.id) == added_ids.end()) {
                                        added.push_back(p);
                                        added_ids.push_back(p.id);
                                    }
                                }
                            }

                            for(Polygon p : added) {
                                cout << "PUSHING: " << p.id << endl;
                                w_polygons.push_back(p);
                            }
                        }
                    }

                    return w_polygons;
                }
            }

            /* add the cross if it didn't cross any w_segments */
            w_segments.push_back(crosses[i]);
        }
    }

    return w_polygons;
}

/* creates a normalized bijection vector to a
 * @param L, the tested segment
 * @param points, the set of datapoints
 * @param n, the size of the points array
 * @return true if the segment is validated, false otherwise
 */
bool test_w_segment(Vector L, double interval, Point *points, int n) {
    int i = 0;
    int j = 0;
    int k = 0;
    Point *t_points = new Point [n]; // holds test points sorted by distance
    vector<double> angles; // holds angles corresponding to test set
    double epsilon = 0.000001;
    Vector V1;
    Vector V2;

    /* find test point set */
    for(k = 0; k < n; k++) {
        t_points[k] = Point(points[k]); // copy data
    }

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

    /* find the midpoint of L */
    Point m = Point(L.start.x, L.start.y, -1);
    m.offset(L.i / 2, L.j / 2);

    /* store radius */
    double radius = distance_p(L.start, m);

    /* run tests on point set */
    for(i = 0; i < n; i++) {
        Point p = Point(points[i]);

        for(j = 0; j < n; j++) {
            Point q = Point(t_points[j]);

            /* skip if points are the same */
            if(p.equals(q)) {
                continue;
            }

            /* skip starting at L.start and L.end */
            if(p.equals(L.start) || p.equals(L.end) || q.equals(L.start) || q.equals(L.end)) {
                continue;
            }

            //printf("CHECKING L (%d, %d) BY %d and %d ON %lf\n", L.start.index, L.end.index, p.index, q.index, interval);

            Vector P = Vector("P", L.start, p);
            Vector Q = Vector("Q", L.start, q);
            Vector P1 = Vector("P1", p, L.start);
            Vector P2 = Vector("P2", p, L.end);
            Vector Q1 = Vector("Q1", q, L.start);
            Vector Q2 = Vector("Q2", q, L.end);

            /* skip if p and q don't straddle L */
            if((determinant(L, P) * determinant(L, Q)) >= 0) {
                //printf("DOES NOT STRADDLE!\n");
                continue;
            }

            /* find intersection between E and L */
            Vector E = Vector("E", p, q);
            Point t = find_intersection(E, L);
            if(t.index == -1) {
                //printf("NO INTERSECTION!\n");
                continue;
            }

            // /* store midpoint of <p, q> */
            // Vector S = Vector("S", p, q);
            // Point s = Point(S.start.x, S.start.y, -1);
            // s.offset(S.i / 2, S.j / 2);
            //
            // /* store grad(f(s)) */
            // double radius = distance_p(m, L.start);
            // Vector G = circular_gradient(m, radius, s);
            // G.normalize();
            // G.i *= radius;
            // G.j *= radius;
            //
            // /* find intersection point between G and circle M */
            // Point g = Point(m);
            // g.offset(G.i, G.j);
            //
            // /* check if points should be tested for ranges */
            // if(abs(distance_p(s, p) - distance_p(s, g)) <= interval) {
            //     //printf("CIRCLE NOT CONTAINED IN M!\n");
            //     continue;
            // }

            /* check if the points are within the circular bijection range */
            if((distance_p(p, m) <= radius) && (distance_p(q, m) <= radius)) {
                //printf("INVALID!\n");
                return false; // invalid
            }
            else {
                //printf("NOT IN BIJECTION RANGE!\n");
            }

         }
    }

    return true; //valid
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

/* generates all shortest paths from W
 * @return paths, the vector of shortest paths
 */
vector<Polygon> generate_final_paths(vector<int *> segments, Point *points, int n) {
    vector<Polygon> paths;
    int *segment = segments[0];

    paths.push_back(dijkstra_polygon(segments, segment, points, n));

    return paths;
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

/* returns the point with the smallest calculated tao-distance
 * @param L, the weslean line being tested
 * @param k, the point generation interval
 * @param points, the array of datapoints
 * @param n, the size of the array
 * @return ..., the line segement to add, either a line or NULL*/
vector<Point> generate_w_points(vector<Point> w_points, Vector L, double interval) {
    /* base case */
    if(L.length < interval) {
        return w_points;
    }
    /* recursive step */
    Point M = Point(L.end.x, L.end.y, -1);
    M.offset(-L.i/2, -L.j/2);
    w_points.push_back(M);
    Vector H1 = Vector("H1", L.start, M);
    w_points = generate_w_points(w_points, H1, interval);
    Vector H2  = Vector("H2", M, L.end);
    w_points = generate_w_points(w_points, H2, interval);
    /* return result */
    return w_points;
}

Point minimum_tao_distance(Vector V, Point *points, int size) {
    double tao;
    double theta;
    double curvature;
    double *tao_dist = new double [size];
    Vector T;
    int i = 0;
    for(i = 0; i < size; i++) {
        Point p = Point(points[i]);
        /* skip start */
        if(p.equals(V.start)) {
            tao_dist[i] = DBL_MAX;
            continue;
        }

        /* initialize vector T */
        T = Vector("T", V.start, p);
        T.normalize();

        /* store tao, theta, and curvature */
        tao = (dot_product(V, T)); //length of V and T is always 1
        if(tao <= -1.0) {
            tao = -1.0;
        }
        else if(tao >= 1.0) {
            tao = 1.0;
        }
        theta = angle_t(tao);
        curvature = calculate_curvature(V, T, tao);
        tao_dist[i] = Vector("T'", V.start, p).length;
    }
    /* find the least tao_distance */
    double curr = tao_dist[0];
    int k = 0;
    for(i = 1; i < size; i++) {
        if(curr > tao_dist[i]) {
            curr = tao_dist[i];
            k = i;
        }
    }
    return points[k];
}

/* removes crossing segments
 * @param segments, the vector of segments
 * @param s, the index of the test segment
 * @param points, the array of points the vector corresponds to
 * @param size, the size of the point array
 * @return the updated segments vector
 */
vector<int *> remove_crossing_segments(vector<int *> segments, int s, Point *points, int size) {
    vector<int *> crosses;
    Vector L = Vector("L", points[segments[s][0]], points[segments[s][1]]);
    int i = 0;
    int j = 0;
    int k = 0;

    /* initialize intersections vector */
    for(i = segments.size() - 1; i >= 0; i--) {
        if(i == s) {
            crosses.push_back(segments[i]);
            continue;
        }
        Vector V = Vector("V", points[segments[i][0]], points[segments[i][1]]);
        if(intersection(V, L)) {
            crosses.push_back(segments[i]);
        }
    }

    /* short circuit */
    if(crosses.size() <= 1) {
        return segments;
    }

    vector<int *> **edge_map = new vector<int *> *[crosses.size()]; // edge_map[i][0,1] = vector<int *>
    double *min_map = new double [crosses.size()];
    vector<int *> edges;
    vector<int> remove;
    int max = 0;
    int index;

    /* store all the edges of crosses[i][0,1], and find the max */
    for(i = 0; i < crosses.size(); i++) {
        edge_map[i] = new vector<int *> [2];

        /* find edge count of the start point */
        edges = edge_search(segments, points[crosses[i][0]].index, points, size);
        for(j = 0; j < crosses.size(); j++) {
            if((index = segment_match(edges, crosses[j][0], crosses[j][1])) > -1) {
                edges.erase(edges.begin() + index);
            }
        }
        edge_map[i][0] = edges;
        edges.clear();
        remove.clear();

        /* find edge count of the end point */
        edges = edge_search(segments, points[crosses[i][1]].index, points, size);
        for(j = 0; j < crosses.size(); j++) {
            if((index = segment_match(edges, crosses[j][0], crosses[j][1])) > -1) {
                edges.erase(edges.begin() + index);
            }
        }
        edge_map[i][1] = edges;
        edges.clear();
        remove.clear();

        min_map[i] = distance_p(points[crosses[i][0]], points[crosses[i][1]]);

        /* store the number of edges and the max */
        if(max < edge_map[i][0].size() + edge_map[i][1].size()) {
            max = edge_map[i][0].size() + edge_map[i][1].size();
        }
    }

    /* find all segments with the max number outgoing edges */
    vector<int> equal;
    for(i = 0; i < crosses.size(); i++) {
        if(max == (edge_map[i][0].size() + edge_map[i][1].size())) {
            equal.push_back(i); // keeps track of which indices are equals
        }
    }

    /* delete all segments that have less edges, and if more than
     * one segment has equal edges delete them all */
    if(equal.size() > 1) {
        /* find the min distance out of the equal segments */
        double min = DBL_MAX;
        for(i = 0; i < equal.size(); i--) {
            if(min > min_map[equal[i]]) {
                min = min_map[equal[i]];
            }
        }
        for(i = 0; i < equal.size(); i--) {
            if(min < min_map[equal[i]]) {
                if((index = segment_match(segments, crosses[equal[i]][0], crosses[equal[i]][1])) > -1) {
                    segments.erase(segments.begin() + index);
                }
            }
        }
    }
    for(i = 0; i < crosses.size(); i++) {
        if(max >= (edge_map[i][0].size() + edge_map[i][1].size())) {
            if((index = segment_match(segments, crosses[i][0], crosses[i][1])) > -1) {
                segments.erase(segments.begin() + index);
            }
        }
    }

    return segments;
}

/* fixes any overlaps with segments
 * @param test, the test segment to check
 * @param segments, the vector of segments to check
 * @param points, the array of datapoins
 * @return the updated segments vector
 */
vector<int *> fix_overlap(int *test, vector<int *> segments, Point *points) {
    int i = 0;
    int j = 0;

    /* construct original segment vector */
    Vector V = Vector("V", points[test[0]], points[test[1]]);

    /* store the segments that overlap and remove them from segments */
    vector<int*> popped;
    for(i = 0; i < segments.size(); i++) {
        Vector T = Vector("T", points[segments[i][0]], points[segments[i][1]]);
        if(overlap(V, T)) {
            popped.push_back(segments[i]);
        }
    }

    popped.push_back(test);

    for(i = 0; i < popped.size(); i++) {
        int index;
        if((index = segment_match(segments, popped[i][0], popped[i][1])) > -1) {
            segments.erase(segments.begin() + index);
        }
    }

    /* create a vector of indices from the popped segments */
    vector<int> sorted;
    for(i = 0; i < popped.size(); i++) {
        sorted.push_back(popped[i][0]);
        sorted.push_back(popped[i][1]);
    }

    /* bubble sort points (indices) by length */
    for(i = 0; i < sorted.size(); i++) {
        for(j = sorted.size() - 1; j > i; j--) {
            double curr;
            double next;
            if(V.i == 0) { // use y values to compare
                curr = points[sorted[j]].y;
                next = points[sorted[j - 1]].y;
            }
            else { // use x values to compare
                curr = points[sorted[j]].x;
                next = points[sorted[j - 1]].x;
            }
            if(curr < next) {
                int tmp = sorted[j];
                sorted[j] = sorted[j - 1];
                sorted[j - 1] = tmp;
            }
        }
    }

    int k;

    /* create the new segments and add them to the vector */
    for(i = 0; i < sorted.size() - 1; i++) {
        if(sorted[i] != sorted[i + 1]) {
            int *tmp_segment = new int[2];
            tmp_segment[0] = sorted[i];
            tmp_segment[1] = sorted[i + 1];
            segments.push_back(tmp_segment);
        }
    }

    return segments;
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

/* creates two polygons, forward and backwards
 * @param edge, the current edge used to create the polygon
 * @param segments, the vector of segments calculated before
 * @param points, the array of datapoints
 * @param size, the size of the point array
 * @return the generated polygons
 */
vector<Polygon> create_polygons(int *input_edge, vector<int *> segments, Point *points, int size) {
    int i = 0;
    int *tmp_segment;
    int *edge = new int [2]; // initialize copy
    edge[0] = input_edge[0];
    edge[1] = input_edge[1];
    int e = edge[0];
    int *tmp_edge = new int [2];
    tmp_edge[0] = edge[0];
    tmp_edge[1] = edge[1];
    vector<int *> processed; // holds processed segments
    vector<Polygon> polygons;

    /* initialize visited array */
    int *visited = new int [size];
    for(i = 0; i < size; i++) {
        visited[i] = 0;
    }

    /* visit the starting point */
    visited[edge[1]] = 1;

    Polygon polygon = Polygon();
    polygon.shape.push_back(edge[0]);

    /* loop to create the polygon
     * 1. using E and V, find the left-most segment
     *    (that's not exactly 180 degrees from Y)
     * 2. set variables and continue with the new segment
     * 3. continue processing until the first point is reached
     */
    vector<int *> edges;
    while(edge[1] != e) { // use e as a test value
        /* add edge(s) to processed segments and polygon*/
        tmp_segment = new int [2];
        tmp_segment[0] = edge[0];
        tmp_segment[1] = edge[1];
        processed.push_back(tmp_segment);
        polygon.shape.push_back(edge[1]);

        /* initialize test variables */
        double max = -180.0;
        int index = -1;
        Vector E = Vector("E", points[edge[0]], points[edge[1]], -1);

        /* calculate edges off of edge[1] */
        edges = edge_search(segments, points[edge[1]].index, points, size);

        /* find the left-most segment */
        for(i = 0; i < edges.size(); i++) {
            if((visited[edges[i][1]] == 1) || (edges[i][1] == edge[0])) { // skip visited points
                continue;
            }
            Vector V = Vector("V", points[edges[i][0]], points[edges[i][1]]);
            //V.normalize();

            /* test if the vector is on the left side */
            if(determinant(E, V) >= 0) {
                /* test if the vector is possibly the greatest */
                if((angle(E, V) * 180 / M_PI) > max && (angle(E, V) * 180 / M_PI) < 180) {
                    max = angle(E, V) * 180 / M_PI;
                    index = edges[i][1];
                }
            }
            /* test if the vector is on the right side */
            else {
                /* test if the vector is possibly the greatest */
                if((angle(E, V) * -180 / M_PI) > max && (angle(E, V) * 180 / M_PI) < 180) {
                    max = angle(E, V) * -180 / M_PI;
                    index = edges[i][1];
                }
            }
        }

        /* pop the point off if the left-most segment was not found */
        if(index == -1) {
            polygon.shape.pop_back();
            edge[1] = edge[0];
            edge[0] = polygon.shape.back();
            /* test if the shape has been popped off to nothing */
            if(edge[0] == edge[1]) {
                polygon = Polygon();
                break;
            }
        }
        else {
            /* set variables */
            visited[index] = 1;
            edge[0] = edge[1];
            edge[1] = index;
        }
    }

    if(polygon.shape.size() >= 2) {
        tmp_segment = new int [2];
        tmp_segment[0] = edge[0];
        tmp_segment[1] = edge[1];
        processed.push_back(tmp_segment);
        polygon.shape.push_back(edge[1]);
        polygon = Polygon(polygon.shape, processed, points); // add all processed segments
        polygons.push_back(polygon);
    }

    edge[0] = tmp_edge[1];
    edge[1] = tmp_edge[0];
    e = edge[0];
    processed.clear();

    /* initialize visited array */
    visited = new int [size];
    for(i = 0; i < size; i++) {
        visited[i] = 0;
    }

    /* visit the starting point */
    visited[edge[1]] = 1;

    polygon = Polygon();
    polygon.shape.push_back(edge[0]);

    /* loop to create the polygon
     * 1. using E and V, find the left-most segment
     *    (that's not exactly 180 degrees from Y)
     * 2. set variables and continue with the new segment
     * 3. continue processing until the first point is reached
     */
    while(edge[1] != e) { // use e as a test value
        /* add edge(s) to processed segments and polygon*/
        tmp_segment = new int [2];
        tmp_segment[0] = edge[0];
        tmp_segment[1] = edge[1];
        processed.push_back(tmp_segment);
        polygon.shape.push_back(edge[1]);

        /* initialize test variables */
        double max = -180.0;
        int index = -1;
        Vector E = Vector("E", points[edge[0]], points[edge[1]], -1);

        /* calculate edges off of edge[1] */
        edges = edge_search(segments, points[edge[1]].index, points, size);

        /* find the left-most segment */
        for(i = 0; i < edges.size(); i++) {
            if((visited[edges[i][1]] == 1) || (edges[i][1] == edge[0])) { // skip visited points
                continue;
            }
            Vector V = Vector("V", points[edges[i][0]], points[edges[i][1]]);
            //V.normalize();

            /* test if the vector is on the left side */
            if(determinant(E, V) >= 0) {
                /* test if the vector is possibly the greatest */
                if((angle(E, V) * 180 / M_PI) > max && (angle(E, V) * 180 / M_PI) < 180) {
                    max = angle(E, V) * 180 / M_PI;
                    index = edges[i][1];
                }
            }
            /* test if the vector is on the right side */
            else {
                /* test if the vector is possibly the greatest */
                if((angle(E, V) * -180 / M_PI) > max && (angle(E, V) * 180 / M_PI) < 180) {
                    max = angle(E, V) * -180 / M_PI;
                    index = edges[i][1];
                }
            }
        }

        /* pop the point off if the left-most segment was not found */
        if(index == -1) {
            polygon.shape.pop_back();
            edge[1] = edge[0];
            edge[0] = polygon.shape.back();
            /* test if the shape has been popped off to nothing */
            if(edge[0] == edge[1]) {
                polygon = Polygon();
                break;
            }
        }
        else {
            /* set variables */
            visited[index] = 1;
            edge[0] = edge[1];
            edge[1] = index;
        }
    }

    if(polygon.shape.size() >= 2) {
        tmp_segment = new int [2];
        tmp_segment[0] = edge[0];
        tmp_segment[1] = edge[1];
        processed.push_back(tmp_segment);
        polygon.shape.push_back(edge[1]);
        polygon = Polygon(polygon.shape, processed, points); // add all processed segments
        polygons.push_back(polygon);
    }

    return polygons;
}

/* searches through a vector of shapes for a duplicate
 * @param polygons, the current vector of polygons
 * @param points, the array of datapoints for printing
 * @return the updated polygon vector
 */
vector<Polygon> delete_duplicate_polygons(vector<Polygon> polygons, Point *points) {
    vector<Polygon> tmp = polygons;
    int i = 0;
    int j = 0;
    int k = 0;

    /* sort each entry */
    for(i = 0; i < tmp.size(); i++) {
        sort((tmp[i]).shape.begin(), (tmp[i]).shape.end());
    }

    /* remove the duplicate from the sorted list */
    for(i = 0; i < tmp.size(); i++) {
        for(j = 0; j < (tmp[i]).shape.size(); j++) {
            for(k = 0; k < (tmp[i]).shape.size(); k++) {
                if(j != k && (tmp[i]).shape[j] == (tmp[i]).shape[k]) {
                    (tmp[i]).shape.erase((tmp[i]).shape.begin() + j);
                    i = 0;
                    j = 0;
                    k = 0;
                }
            }
        }
    }

    /* check for an equal entry */
    for(i = 0; i < tmp.size(); i++) {
        vector<int> remove;
        for(j = tmp.size(); j >= 0; j--) {
            if(((tmp[i]).shape.size() == (tmp[j]).shape.size()) && (i != j)) {
                bool complete = true;
                for(k = 0; k < (tmp[i]).shape.size(); k++) {
                    if((tmp[i]).shape[k] != (tmp[j]).shape[k]) {
                        complete = false;
                        break;
                    }
                }
                if(complete) {
                    if (find(remove.begin(), remove.end(), j) == remove.end()) {
                        remove.push_back(j);
                    }
                }
            }
        }
        /* delete entries */
        for(j = 0; j < remove.size(); j++) {
            polygons.erase(polygons.begin() + remove[j]);
            tmp.erase(tmp.begin() + remove[j]);
            i = 0;
        }
    }

    return polygons;
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
