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

    /* check if the segment is found */
    for(i = 0; i < segments.size(); i++) {
        /* find the indexed value */
        for(j = 0; j < size; j++) {
            if(points[j].index == vertex) {
                break;
            }
        }
        if(segments[i][0] == points[vertex].index) {
            tmp = new int [2];
            tmp[0] = j;
            for(j = 0; j < size; j++) {
                if(points[j].index == points[segments[i][1]].index) {
                    break;
                }
            }
            tmp[1] = j;
            edges.push_back(tmp);
        }
        else if(segments[i][1] == points[vertex].index) {
            tmp = new int [2];
            tmp[0] = j;
            for(j = 0; j < size; j++) {
                if(points[j].index == points[segments[i][0]].index) {
                    break;
                }
            }
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

/* returns true if the points intersect, false otherwise */
bool intersection(Vector V1, Vector V2)
{
    Point p1 = V1.start;
    Point p2 = V1.end;
    Point p3 = V2.start;
    Point p4 = V2.end;
    Point tmp;
    double y = 0.0;
    double x = 0.0;
    double m1 = 0.0;
    double m2 = 0.0;
    double b1 = 0.0;
    double b2 = 0.0;
    double r = DBL_MAX;
    double l = DBL_MIN;
    /* find intersection area */
    if (p2.x < p1.x) {
        tmp = p1;
        p1 = p2;
        p2 = tmp;
    }
    if (p4.x < p3.x) {
        tmp = p3;
        p3 = p4;
        p4 = tmp;
    }
    /* make sure none of the nodes are the same */
    if((p1.index == p3.index) || (p1.index == p4.index) || (p2.index == p3.index) || (p2.index == p4.index)) {
        return false;
    }
    /* find the left and right boundaries of intersection on the x-axis */
    if(p1.x >= p3.x) {
        /*           p1       */
        /*            .       */
        /*            .       */
        /*      p3---------p4 */
        /*            .       */
        /*            .       */
        /*            p2      */
        if(p2.x == p1.x) {
            p1.x -= 0.000001;
            p2.x += 0.000001;
            p3.x -= 0.000001;
            p4.x += 0.000001;
        }
        l = p1.x;
        /* p3--------------------p4 */
        /*      p1---------p2 */
        if(p2.x <= p4.x) {
            r = p2.x;
        }
        /* p3-----------p4 */
        /*        p1-----------p2 */
        else {
            r = p4.x;
        }
    }
    else if(p2.x >= p3.x) {
        /*           p3       */
        /*            .       */
        /*            .       */
        /*      p1---------p2 */
        /*            .       */
        /*            .       */
        /*            p4      */
        if(p3.x == p4.x) {
            p1.x -= 0.000001;
            p2.x += 0.000001;
            p3.x -= 0.000001;
            p4.x += 0.000001;
        }
        l = p3.x;
        /* p1--------------------p2 */
        /*      p3---------p4 */
        if(p2.x >= p4.x) {
            r = p4.x;
        }
        /* p1-----------p2 */
        /*        p3-----------p4 */
        else {
            r = p2.x;
        }
    }
    /* p1---------p2  p3---------p4 */
    else {
        return false;
    }
    /* p3---------p4  p1---------p2 */
    if(l >= r) {
        return false;
    }
    /* create first equation */
    y = p1.y;
    m1 = (p2.y - p1.y) / (p2.x - p1.x);
    x = p1.x;
    b1 = y - m1 * x;
    /* create second equation */
    y = p3.y;
    m2 = (p4.y - p3.y) / (p4.x - p3.x);
    x = p3.x;
    b2 = y - m2 * x;
    /* solve linear system */
    x = (b2 - b1) / (m1 - m2);
    y = m1 * x + b1;
    /* don't push the segment if it is an overlap */
    if((b2 == b1) && (m1 == m2)) {
        if(p3.x < p2.x) {
            return true;
        }
        else if(p1.x < p4.x) {
            return true;
        }
    }
    /* don't push the segment if there is an intersection */
    else if(x <= r && x >= l) {
        return true;
    }
    return 0;
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

/* calculates distance given index and structure */
double tao_distance(Vector V, double curvature, double theta)
{
    return (V.length + curvature + theta);
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

/* prints curvature structure for debugging */
void print(Vector V, Vector T1, Vector T2, double curvature, double theta, double tao, double tao_distance)
{
    printf("V: %d[0](%lf, %lf), %d[1](%lf, %lf), <%lf, %lf>, |V| = %lf\n", V.start.index, V.start.x, V.start.y, V.end.index, V.end.x, V.end.y, V.i, V.j, V.length);
    printf("T1: point[0](%lf, %lf), point[1](%lf, %lf), <%lf, %lf>, |T1| = %lf\n", T1.start.x, T1.start.y, T1.end.x, T1.end.y, T1.i, T1.j, T1.length);
    printf("T2: point[0](%lf, %lf), point[1](%lf, %lf), <%lf, %lf>, |T2| = %lf\n", T2.start.x, T2.start.y, T2.end.x, T2.end.y, T2.i, T2.j, T2.length);
    printf("curvature: %lf; ", curvature);
    printf("angle = %lf; ", theta * 180 / M_PI);
    printf("tao = %lf; ", tao);
    printf("tao-distance = %lf\n\n", tao_distance);
}

/* prints to the terminal if there is an error assigning memory */
void memory_error(void)
{
    printf("\n\nError assigning memory. Exiting Program. Good Day.\n\n");
}

/* experimental algorithm... */
vector<Polygon> w_polygon_construction(Point *points, int size, FILE *gnu_files[NUM_FILES]) {

    /////////////////////////////////
    // INITIALIZATION OF VARIABLES //
    /////////////////////////////////

    vector<int *> segments;
    int *tmp_segment = new int [2];
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;

    ////////////////////////
    // START CALCULATIONS //
    ////////////////////////

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

    interval /= 2;
    vector<int *> crosses;

    /* test all line segments */
    for(i = 0; i < size; i++) {
        for(j = 0; j < size; j++) {
            if(i == j) {
                continue;
            }
            Vector L = Vector("L", points[i], points[j]);
            if(segment_match(segments, i, j) == -1) {
                /* if the line is valid, add it as a segment */
                if(test_w_segment(L, interval, points, size)) {
                    tmp_segment = new int [2];
                    tmp_segment[0] = i;
                    tmp_segment[1] = j;

                    /* record segment */
                    if(segments.size() == 0) {
                        segments.push_back(tmp_segment);
                    }
                    else {
                        segments = fix_overlap(tmp_segment, segments, points);
                    }
                }
            }
        }
    }

    /* bubble sort segments by length, to aid in removing correct crosses */
    for(i = 0; i < segments.size(); i++) {
        for(j = segments.size() - 1; j > i; j--) {
            double curr = distance_p(points[segments[j][0]], points[segments[j][1]]);
            double next = distance_p(points[segments[j - 1][0]], points[segments[j - 1][1]]);
            if(curr < next) {
                int *tmp = new int [2];
                tmp[0] = segments[j][0];
                tmp[1] = segments[j][1];
                segments[j][0] = segments[j - 1][0];
                segments[j][1] = segments[j - 1][1];
                segments[j - 1][0] = tmp[0];
                segments[j - 1][1] = tmp[1];
            }
        }
    }

    /* remove crossing segments */
    int prev_size = segments.size();
    i = 0;
    while(i < segments.size()) {
        segments = remove_crossing_segments(segments, i, points, size);
        if(prev_size == segments.size()) {
            i++;
        }
        else {
            i = 0;
            prev_size = segments.size();
        }
    }

    printf("\nFILTERED SEGMENTS:\n");
    for(i = 0; i < segments.size(); i++) {
        printf("%d: (%d, %d)\n", i, points[segments[i][0]].index, points[segments[i][1]].index);
    }

    /* create w_polygons */
    vector<Polygon> polygons;
    vector<int *> edges;
    for(i = 0; i < size; i ++) {
        /* find the polygon starting at each edge */
        edges = edge_search(segments, points[i].index, points, size);
        for(int *edge : edges) {
            Polygon tmp_polygon = create_polygon(edge, segments, points, size);
            /* push the polygon */
            if(tmp_polygon.shape.size() > 0) {
                polygons.push_back(tmp_polygon);
            }
        }
    }

    return polygons;
}

/* tests if a weslean segment is valid
 * @param L, the weslean line being tested
 * @param interval, the weslean point generation interval
 * @param points, the array of datapoints
 * @param n, the size of the array
 * @return ..., the line segement to add, either a line or NULL*/
bool test_w_segment(Vector L, double interval, Point *points, int n) {
    bool valid = true;
    int i = 0;
    int j = 0;
    int k = 0;

    //printf("TESTING: (%d, %d)\n", L.start.index, L.end.index);

    /* create initial point set */
    int l = 0;
    int m = 0;
    Point *t_points = new Point [n + 1];
    for(; l < n; l++) {
        t_points[m++] = points[l];
    }

    /* create w_points */
    vector<Point> w_points;
    w_points = generate_w_points(w_points, L, interval);

    /* test all w_points, q = {q_1, q_2, ..., q*} */
    for(j = 0; j < w_points.size(); j++) {
        Point q = Point(w_points[j]);
        bool skip = false;
        /* skip q if it is equal to any of the actual points */
        for(i = 0; i < n; i++) {
            if(q.equals(points[i])) {
                skip = true;
                break;
            }
        }
        if(skip) {
            continue;
        }
        /* add q to the point set */
        t_points[m] = q;

        /* test point set */
        for(i = 0; i < n + 1; i++) {
            /* skip starting at q, L.start, and L.end */
            if(t_points[i].equals(q) || t_points[i].equals(L.start) || t_points[i].equals(L.end)) {
                continue;
            }

            Point p = Point(t_points[i]);

            /* find the forward connection, i.e.
             * find the minimum tao_distance between all
             * vectors <p, *> and Q = <p, q>
             */
            Vector Q = Vector("Q", p, q);
            Q.normalize();
            Point t = minimum_tao_distance(Q, t_points, n + 1);
            //printf("%d: %d\n", p.index, t.index);

            /* check if the line needs further testing */
            if(!t.equals(q)) {
                continue;
            }

            /* find the follow connection, i.e.
             * find the minimum tao_distance between all
             * vectors <q, p> and Q = <q, q'>
             */
            Q = Vector("Q", p, q, -1); // shift and normalize
            l = 0;
            Point r;
            if(L.start.y == L.end.y) {
                if(L.start.x <= L.end.x) {
                    r = Point(L.start);
                }
                else {
                    r = Point(L.end);
                }
            }
            else if(L.start.y < L.end.y) {
                r = Point(L.start);
            }
            else {
                r = Point(L.end);
            }
            Vector R = Vector("R", r, q, -1); // shift and normalize
            //R.print();
            for(k = 0; k < n; k++) {
                Point s = Point(points[k]);
                Vector S = Vector("S", q, s);
                S.normalize();
                /* check if point s should remain in the test set */
                if((determinant(Q, R) < 0 && determinant(S, R) < 0) || (determinant(Q, R) > 0 && determinant(S, R) > 0)) {
                    //S.print();
                    //printf("  determinant = %lf\n", determinant(R, S));
                    //printf("... %d is GOOD!\n", s.index);
                    t_points[l++] = s;
                }
            }
            t_points[l++] = L.start;
            t_points[l++] = L.end;
            //printf("q':<(%0.3lf, %0.3lf), (%0.3lf, %0.3lf)>\n", Q.start.x, Q.start.y, Q.end.x, Q.end.y);
            t = minimum_tao_distance(Q, t_points, l);
            //printf("%d: %d\n", p.index, t.index);

            /* check if the line was found to be invalid */
            if(!t.equals(L.start) && !t.equals(L.end)) {
                valid = false;
                break;
            }
        }
        if(!valid) {
            break;
        }
    }
    return valid;
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
        //tao_dist[i] = tao_distance(Vector("T'", V.start, p), curvature, theta);
        tao_dist[i] = Vector("T'", V.start, p).length;
        //printf("%d -> %d: %0.3lf + %0.3lf + %0.3lf = %0.3lf\n", V.start.index, p.index, distance_p(V.start, p), curvature, theta, tao_dist[i]);
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

    /* construct original segment points and vector */
    Point p1 = Point(points[test[0]]);
    Point p2 = Point(points[test[1]]);
    Vector V = Vector("V", p1, p2);
    if(V.i == 0) { // use y values
        if(p2.y < p1.y) {
            Point tmp = p2;
            p2 = p1;
            p1 = tmp;
            Vector V = Vector("V", p1, p2);
        }
    }
    else { // use x values
        if(p2.x < p1.x) {
            Point tmp = p2;
            p2 = p1;
            p1 = tmp;
            Vector V = Vector("V", p1, p2);
        }
    }

    /* store the segments that overlap and remove them from segments */
    vector<int*> popped;
    for(i = segments.size() - 1; i >= 0; i--) {

        /* construct test segment points and vector */
        Point p3 = Point(points[segments[i][0]]);
        Point p4 = Point(points[segments[i][1]]);
        Vector T = Vector("T", p3, p4);
        if(T.i == 0) { // use y values
            if(p4.y < p3.y) {
                Point tmp = p4;
                p4 = p3;
                p3 = tmp;
                Vector T = Vector("T", p3, p4);
            }
        }
        else { // use x values
            if(p4.x < p3.x) {
                Point tmp = p4;
                p4 = p3;
                p3 = tmp;
                Vector T = Vector("T", p3, p4);
            }
        }

        Point u = Point(p1);
        u.offset(T.i, T.j);
        Vector U = Vector("U", p1, u);
        /* parallel test */
        if(determinant(V, U) == 0) {
            double y = 0.0;
            double x = 0.0;
            double m1 = 0.0;
            double m2 = 0.0;
            double b1 = 0.0;
            double b2 = 0.0;
            /* find which intercept to use, x or y */
            if(V.i == 0) { // use x intercept
                /* create first equation */
                y = p1.x;
                m1 = (p2.x - p1.x) / (p2.y - p1.y);
                x = p1.y;
                b1 = y - m1 * x;
                /* create second equation */
                y = p3.x;
                m2 = (p4.x - p3.x) / (p4.y - p3.y);
                x = p3.y;
                b2 = y - m2 * x;
            }
            else { // use y intercept
                /* create first equation */
                y = p1.y;
                m1 = (p2.y - p1.y) / (p2.x - p1.x);
                x = p1.x;
                b1 = y - m1 * x;
                /* create second equation */
                y = p3.y;
                m2 = (p4.y - p3.y) / (p4.x - p3.x);
                x = p3.x;
                b2 = y - m2 * x;
            }
            if(b1 == b2) {
                //popped.push_back(segments[i]);
                if(V.i == 0) { // use y values
                    if(p1.y <= p3.y) { // P1 < P3
                        if(p3.y <= p2.y) { // P1 --> P3 --> P2 --> P4
                            popped.push_back(segments[i]);
                        }
                    }
                }
                else { // use x values
                    if(p1.x <= p3.x) { // P1 < P3
                        if(p3.x <= p2.x) { // P1 --> P3 --> P2 --> P4
                            popped.push_back(segments[i]);
                        }
                    }
                }
            }
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

/* creates a polygon
 * @param edge, the current edge used to create the polygon
 * @param segments, the vector of segments calculated before
 * @param points, the array of datapoints
 * @param size, the size of the point array
 * @return the generated polygon
 */
Polygon create_polygon(int *edge, vector<int *> segments, Point *points, int size) {
    int i = 0;
    int e = edge[0];
    vector<int> shape;
    Polygon polygon = Polygon();
    polygon.shape.push_back(edge[0]);
    polygon.perimeter = distance_p(points[edge[0]], points[edge[1]]);

    /* loop to create the polygon
     * 1. using E and V, find the left-most segment
     *    (that's not exactly 180 degrees from Y)
     * 2. set variables and continue with the new segment
     * 3. continue processing until the first point is reached
     */
    vector<int *> edges;
    while(edge[1] != e) { // use e as a test value
        polygon.shape.push_back(edge[1]);
        /* initialize test variables */
        double max = 0.0;
        int index = -1;
        Vector E = Vector("E", points[edge[0]], points[edge[1]], -1);
        edges = edge_search(segments, points[edge[1]].index, points, size);

        /* find the left-most segment */
        for(i = 0; i < edges.size(); i++) {
            if(edges[i][1] == edge[0]) { // skip the origin point
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
        }

        /* return null if the left-most segment was not found */
        if(max <= 0.0) {
            Polygon empty;
            return empty;
        }

        /* set variables */
        edge[0] = edge[1];
        edge[1] = index;
        polygon.perimeter += distance_p(points[edge[0]], points[edge[1]]);
    }
    polygon.shape.push_back(edge[1]);
    polygon = Polygon(polygon.shape, points);

    return polygon;
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
