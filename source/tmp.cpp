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

/* runs tao-distance algorithm on dataset and generates optimal segments */
void construct_segments(vector<int *> *segments, Point *points, Point begin, int n, int size, FILE *gnu_files[NUM_FILES], int *mapped, int **recorded)
{
    /////////////////////////////////
    // INITIALIZATION OF VARIABLES //
    /////////////////////////////////

    Vector V;
    Vector T1;
    Vector T2;
    Point *curr = new Point [size];
    Point best;
    Point start;
    Point prev;
    Point center;
    double sum_x = 0.0;
    double sum_y = 0.0;
    int **tmp_segments = new int * [size + 1];
    int *pushed_segment;
    int *loop = new int [size];
    int *visited = new int [size];
    int total_size = size;
    int visited_count = 0;
    int count = 0;
    int i = 0;
    int j = 0;
    int k = 0;
    int m = 0;
    /* initialization of point traversal and segment arrays */
    for(i = 0; i < size; i++) {
        curr[i].x = DBL_MAX;
        curr[i].y = DBL_MAX;
        curr[i].tao_distance = DBL_MAX;
        curr[i].index = -1;
        tmp_segments[i] = new int [2];
        tmp_segments[i][0] = INT_MAX;
        tmp_segments[i][1] = INT_MAX;
    }
    tmp_segments[i] = new int [2];
    tmp_segments[i][0] = INT_MAX;
    tmp_segments[i][1] = INT_MAX;
    best.x = DBL_MIN;
    best.y = DBL_MIN;
    best.tao_distance = DBL_MIN;
    best.index = -1;
    /* initialize center point */
    for(i = 0; i < size; i++) {
        sum_x += points[i].x;
        sum_y += points[i].y;
    }
    center.x = sum_x / size;
    center.y = sum_y / size;

    /* initialize start, prev, and begin */
    start = Point(begin);
    prev = Point(begin);
    best = Point(begin);
    loop[m++] = n;

    /* initialize vector T1 */
    T1 = Vector("T1", center, start, INT_MAX);
    T1.offset(start.x - center.x, start.y - center.y);

    ///////////////////////
    // PLOT CENTER POINT //
    ///////////////////////

    fprintf(gnu_files[4], "%lf %lf\n", center.x, center.y);

    ////////////////////////
    // START CALCULATIONS //
    ////////////////////////

    /* outer loop, calculates total distance
     * condition: loop until all points are visited
     * T'(n):     O(n)
     */
    while(visited_count <= total_size) {
        /* store start index in visited-array */
        visited[start.index] = 1;
        /* refresh loop variables */
        i = 0;
        best.tao_distance = DBL_MAX;
        best.index = start.index;
        /* initialize vector T2 = 0 */
        T2 = Vector("T2", start, start);
        /* initialize vector V = 0 */
        V = Vector("V", start, start);
        count = 0;
        /* loops through all possible indices from start
         * condition:  loop until all points are visited
         * T'(n):      O(n - 1)
         */
        while(count < size) {
            /* skip current index and previous index */
            if((points[i].equals(best)) || (points[i].equals(prev))) {
                curr[i].tao_distance = DBL_MAX;
                i++;
                count++;
                continue;
            }
            /* initialize vector V */
            V.end = Point(points[i]);
            V.refresh();
            /* initialize vector T2 */
            T2 = Vector("T2", T1.start, V.end, INT_MAX);
            /* store tao, theta, and curvature */
            curr[i].tao = (dot_product(T1, T2)); //length of T1 and T2 is always 1
            if(curr[i].tao <= -1.0) {
                curr[i].tao = -1.0;
            }
            else if(curr[i].tao >= 1.0) {
                curr[i].tao = 1.0;
            }
            curr[i].x = V.end.x;
            curr[i].y = V.end.y;
            curr[i].index = V.end.index;
            curr[i].theta = angle_t(curr[i].tao);
            curr[i].curvature = calculate_curvature(T1, T2, curr[i].tao);
            //curr[i].tao_distance = tao_distance(V, curr[i].curvature, curr[i].theta) * curr[i].theta; //addition of multiply by theta
            curr[i].tao_distance = tao_distance(V, curr[i].curvature, curr[i].theta); //addition of multiply by theta
            V.end = Point(curr[i]);
            i++;
            count++;
        }
        /* sets the previous point as the previous best point */
        prev = best;
        /* find point with the lowest tao-distance */
        for(i = 0; i < size; i++) {
            if(best.tao_distance > curr[i].tao_distance) {
                best = Point(curr[i]);
                k = i;
                
            }
        }
        /* displays points of equal tao-distance */
        for(i = 0; i < size; i++) {
            if(best.index == curr[i].index) {
                ;
            }
            else if(best.tao_distance == curr[i].tao_distance) {
                //printf("(%d) %lf = (%d) %lf\n", best.index, best.tao_distance, curr[i].index, curr[i].tao_distance);
            }
        }
        /* record path */
        loop[m++] = k;
        /* if the best point has been visited before */
        if(visited[best.index] == 1) {
            m--;
            if(m != 0) {
                for(j = 0; j < m; j++) {
                    //printf("%d->", points[loop[j]].index);
                    mapped[loop[j]] = 1;
                }
                //printf("%d\n", points[loop[m]].index);
                for(j = 0; j < m; j++) {
                    /* record segment */
                    tmp_segments[j][0] = loop[j];
                    tmp_segments[j][1] = loop[j + 1];
                }
                /* calculates the tmp_segments for each contour */
                for(j = 0; j < m; j++) {
                    /* skips over recorded tmp_segments */
                    if((recorded[tmp_segments[j][0]][tmp_segments[j][1]] == 1) || (recorded[tmp_segments[j][1]][tmp_segments[j][0]] == 1)) {
                        continue;
                    }
                    /* pushes new segments */
                    pushed_segment = new int [2];
                    pushed_segment[0] = tmp_segments[j][0];
                    pushed_segment[1] = tmp_segments[j][1];
                    segments->push_back(pushed_segment);
                    recorded[tmp_segments[j][0]][tmp_segments[j][1]] = 1;
                    /* book-keeping */
                    //printf("%d = (%d, %d): <%d,%d>\n", j, tmp_segments[j][0], tmp_segments[j][1], points[tmp_segments[j][0]].index, points[tmp_segments[j][1]].index);
                }
            }
            return;
        }
        /* reinitializing vector T1 */
        T2 = Vector("T2", T1.start, best, INT_MAX);
        T1 = Vector("T1", best, best);
        T1.end.offset(T2.i, T2.j);
        T1.refresh();
        /* shifts starting point to best point */
        start = Point(best);
        /* initializing vector T2 */
        T2 = Vector("T2", start, start);
        /* initializing vector V */
        V = Vector("V", start, start);
        count = 0;
        visited_count++;
    }
    return;
}

/* calculates the connections for un-joined vertices */
void join_vertex(vector<int *> *segments, Point *points, Point begin, int n, int size)
{
    Vector V;
    Vector T1;
    Vector T2;
    Point *curr = new Point [size];
    Point best;
    Point start;
    Point prev;
    Point center;
    double sum_x = 0.0;
    double sum_y = 0.0;
    int **tmp_segments = new int * [size + 1];
    int *pushed_segment;
    int *visited = new int [size];
    int count = 0;
    int added = 0;
    int i = 0;
    int j = 0;
    int k = 0;
    int m = 0;
    /* initialization of arrays */
    for(i = 0; i < size; i++) {
        curr[i].x = DBL_MAX;
        curr[i].y = DBL_MAX;
        curr[i].tao_distance = DBL_MAX;
        curr[i].index = INT_MAX;
        tmp_segments[i] = new int [2];
        tmp_segments[i][0] = INT_MAX;
        tmp_segments[i][1] = INT_MAX;
    }
    tmp_segments[i] = new int [2];
    tmp_segments[i][0] = INT_MAX;
    tmp_segments[i][1] = INT_MAX;
    best.x = DBL_MAX;
    best.y = DBL_MAX;
    best.tao_distance = DBL_MAX;
    best.index = -1;

    /* calculate average point */
    for(i = 0; i < size; i++) {
        sum_x += points[i].x;
        sum_y += points[i].y;
    }
    center.x = sum_x / size;
    center.y = sum_y / size;

    start = Point(begin);
    prev = Point(begin);
    best = Point(begin);
    /* initializing vector T1 */
    T1 = Vector("T1", center, start, INT_MAX);
    T1.offset(start.x - center.x, start.y - center.y);
    /* outer loop */
    while(m < 2) {
        /* store start index in visited-array */
        visited[start.index] = 1;
        i = 0;
        /* refreshing best index */
        best.tao_distance = DBL_MAX;
        best.index = start.index;
        /* initializing vector T2 */
        T2 = Vector("T2", start, start);
        /* initializing vector V */
        V = Vector("V", start, start);
        count = 0;
        /* loops through all possible indices from start */
        while(count < size) {
            /* skip current index and previous index */
            if((points[i].index == best.index) || (points[i].index == prev.index)) {
                curr[i].tao_distance = DBL_MAX;
                i++;
                count++;
                continue;
            }
            /* initializing vector V */
            V.end = Point(points[i]);
            V.refresh();
            /* initializing vector T2 */
            T2 = Vector("T2", T1.start, V.end, INT_MAX);
            /* initializing tao, theta, and curvature */
            curr[i].tao = (dot_product(T1, T2)); //length of T1 and T2 is always 1
            if(curr[i].tao <= -1.0) {
                curr[i].tao = -1.0;
            }
            else if(curr[i].tao >= 1.0) {
                curr[i].tao = 1.0;
            }
            curr[i].x = V.end.x;
            curr[i].y = V.end.y;
            curr[i].index = V.end.index;
            curr[i].theta = angle_t(curr[i].tao);
            curr[i].curvature = calculate_curvature(T1, T2, curr[i].tao);
            curr[i].tao_distance = tao_distance(V, curr[i].curvature, curr[i].theta) * curr[i].theta; //addition of multiply by theta
            V.end.tao_distance = curr[i].tao_distance;
            i++;
            count++;
        }
        /* sets the previous point as the previous best point */
        prev = best;
        /* find point with the lowest tao-distance */
        for(i = 0; i < size; i++) {
            if(best.tao_distance > curr[i].tao_distance) {
                best = Point(curr[i]);
                k = i;
            }
        }
        /* record segment */
        tmp_segments[m][0] = n;
        tmp_segments[m][1] = k;
        /* pushes new segments */
        pushed_segment = new int [2];
        pushed_segment[0] = tmp_segments[m][0];
        pushed_segment[1] = tmp_segments[m][1];
        segments->push_back(pushed_segment);
        m++;
        /* reinitializing vector T1 */
        T2 = Vector("T2", T1.start, best, INT_MAX);
        T1 = Vector("T1", best, best);
        T1.end.offset(T2.i, T2.j);
        T1.refresh();
        /* shifts starting point to best point */
        start = Point(best);
        /* initializing vector T2 */
        T2 = Vector("T2", start, start);
        /* initializing vector V */
        V = Vector("V", start, start);
        count = 0;
    }
    return;
}

/* calculates the connections for un-joined segments */
void join_segment(vector<int *> *segments, Point *points, Point begin, Point end, int n, int m, int size)
{
    Vector V;
    Vector T1;
    Vector T2;
    Point *curr = new Point [size];
    Point best;
    Point start;
    Point prev;
    Point center;
    double sum_x = 0.0;
    double sum_y = 0.0;
    int **tmp_segments = new int * [size + 1];
    int *pushed_segment;
    int *visited = new int [size];
    int count = 0;
    int added = 0;
    int i = 0;
    int j = 0;
    int k = 0;
    /* initialization of arrays */
    for(i = 0; i < size; i++) {
        curr[i].x = DBL_MAX;
        curr[i].y = DBL_MAX;
        curr[i].tao_distance = DBL_MAX;
        curr[i].index = INT_MAX;
        tmp_segments[i] = new int [2];
        tmp_segments[i][0] = INT_MAX;
        tmp_segments[i][1] = INT_MAX;
    }
    tmp_segments[i] = new int [2];
    tmp_segments[i][0] = INT_MAX;
    tmp_segments[i][1] = INT_MAX;
    best = Point(begin);
    best.tao_distance = DBL_MAX;
    /* initializing vector T1 */
    T1 = Vector("T1", end, begin, INT_MAX);
    T1.offset(begin.x - end.x, begin.y - end.y);
    /* store start index in visited-array */
    visited[begin.index] = 1;
    i = 0;
    /* refreshing best index */
    best.tao_distance = DBL_MAX;
    best.index = begin.index;
    /* initializing vector T2 */
    T2 = Vector("T2", start, start);
    /* initializing vector V */
    V = Vector("V", start, start);
    count = 0;
    /* loops through all possible indices from start */
    while(count < size) {
        /* skip current index */
        if((points[i].index == points[n].index) || (points[i].index == points[m].index)) {
            curr[i].tao_distance = DBL_MAX;
            i++;
            count++;
            continue;
        }
        /* initializing vector V */
        V.end = Point(points[i]);
        V.refresh();
        /* initializing vector T2 */
        T2 = Vector("T2", T1.start, V.end, INT_MAX);
        /* initializing tao, theta, and curvature */
        curr[i].tao = (dot_product(T1, T2)); //length of T1 and T2 is always 1
        if(curr[i].tao <= -1.0) {
            curr[i].tao = -1.0;
        }
        else if(curr[i].tao >= 1.0) {
            curr[i].tao = 1.0;
        }
        curr[i].x = V.end.x;
        curr[i].y = V.end.y;
        curr[i].index = V.end.index;
        curr[i].theta = angle_t(curr[i].tao);
        curr[i].curvature = calculate_curvature(T1, T2, curr[i].tao);
        curr[i].tao_distance = tao_distance(V, curr[i].curvature, curr[i].theta) * curr[i].theta; //addition of multiply by theta
        V.end.tao_distance = curr[i].tao_distance;
        i++;
        count++;
    }
    /* find point with the lowest tao-distance */
    for(i = 0; i < size; i++) {
        if(best.tao_distance > curr[i].tao_distance) {
            best = Point(curr[i]);
            k = i;
        }
    }
    /* record segment */
    tmp_segments[m][0] = m;
    tmp_segments[m][1] = k;
    //printf("... joining: <%d,%d>\n", points[n].index, points[k].index);
    /* pushes new segments */
    pushed_segment = new int [2];
    pushed_segment[0] = tmp_segments[m][0];
    pushed_segment[1] = tmp_segments[m][1];
    segments->push_back(pushed_segment);
    return;
}

vector<struct polygon_t> construct_polygons(vector<int *> segments, Point *points, int size)
{
    vector<vector<int> > tessellations;
    vector<struct polygon_t> polygons;
    struct polygon_t polygon;
    double sum_x = 0.0;
    double sum_y = 0.0;
    int i = 0;
    /* traverse the i-th branch of each point */
    for(i = 0; i < size; i++) {
        /* right-right addition */
        tessellations = tessellate(tessellations, segments, points, size, 'r', 'r', i);
        /* right-left addition */
        tessellations = tessellate(tessellations, segments, points, size, 'r', 'l', i);
        /* left-left addition */
        tessellations = tessellate(tessellations, segments, points, size, 'l', 'l', i);
        /* left-right addition */
        tessellations = tessellate(tessellations, segments, points, size, 'l', 'r', i);
    }
    /* stores in polygon_t structure format */
    for(i = 0; i < tessellations.size(); i++) {
        polygon.shape = tessellations[i];
        polygon.perimeter = find_perimeter(tessellations[i], points);
        polygons.push_back(polygon);
    }
    return polygons;
}

/* calculate tessellations of polygons given all contoured segments */
vector<vector<int> > tessellate(vector<vector<int> > tessellations, vector<int *> segments, Point *points, int size, char init, char add, int branch)
{
    vector<int *> edges; //pointer contains the index followed by position
    vector<int> *shape = NULL;
    Vector X; //X-axis vector
    Vector Y; //Y-axis vector
    Point start;
    int i = 0;

    for(i = 0; i < size; i++) {
        /* find the initial cluster of edges */
        edges = edge_search(segments, points[i].index, points, size);
        /* skip if the requested branch is out of bounds */
        if(edges.size() < branch) {
            break;
        }
        start = points[i];
        Y = Vector("Y", start, start);
        Y.end.offset(0, 1);
        Y.refresh();
        X = Vector("X", start, start);
        X.end.offset(1, 0);
        X.refresh();
        shape = find_shape(segments, points, start, size, init, add, X, Y);
        /* skip NULL shapes */
        if(!shape) {
            continue;
        }
        tessellations.push_back(*shape);
    }
    return tessellations;
}

/* finds the four shapes of a cross, given the indices of the crossing segments */
vector<struct polygon_t> tessellate_cross(vector<int *> segments, int i, int j, Point *points, int size)
{
    vector<struct polygon_t> tessellations;
    struct polygon_t polygon;
    Vector X; //X-axis vector
    Vector Y; //Y-axis vector
    Point start;
    Point prev;
    vector<int> *shape = NULL;
    vector<int *> segments_i = segments;
    segments_i.erase(segments_i.begin() + j);
    vector<int *> segments_j = segments;
    segments_j.erase(segments_j.begin() + i);
    if(i > j) {
        i--;
    }
    else {
        j--;
    }
    /*
    int k = 0;
    printf("I-SEGMENTS:\n");
    for(k = 0; k < segments_i.size(); k++) {
        printf("(%d, %d)\n", points[segments_i[k][0]].index, points[segments_i[k][1]].index);
    }
    printf("J-SEGMENTS:\n");
    for(k = 0; k < segments_j.size(); k++) {
        printf("(%d, %d)\n", points[segments_j[k][0]].index, points[segments_j[k][1]].index);
    }*/
    /* find the first shape using segments[i][0] */
    start = Point(points[segments_i[i][0]]);
    prev = Point(points[segments_i[i][1]]);
    Y = Vector("Y", prev, start, INT_MAX);
    Y.offset(start.x - prev.x, start.y - prev.y);
    X = Vector("X", start, start, INT_MAX);
    X.end.offset(Y.j, -Y.i);
    X.refresh();
    shape = find_shape(segments_i, points, start, size, 'r', 'r', X, Y);
    /* skip NULL shapes */
    if(shape) {
        polygon.shape = *shape;
        polygon.perimeter = find_perimeter(*shape, points);
        tessellations.push_back(polygon);
    }
    /* find the second shape using segments[i][1] */
    start = points[segments_i[i][1]];
    prev = points[segments_i[i][0]];
    Y = Vector("Y", prev, start, INT_MAX);
    Y.offset(start.x - prev.x, start.y - prev.y);
    X = Vector("X", start, start, INT_MAX);
    X.end.offset(Y.j, -Y.i);
    X.refresh();
    shape = find_shape(segments_i, points, start, size, 'r', 'r', X, Y);
    /* skip NULL shapes */
    if(shape) {
        polygon.shape = *shape;
        polygon.perimeter = find_perimeter(*shape, points);
        tessellations.push_back(polygon);
    }
    /* find the third shape using segments[j][0] */
    start = points[segments_j[j][0]];
    prev = points[segments_j[j][1]];
    Y = Vector("Y", prev, start, INT_MAX);
    Y.offset(start.x - prev.x, start.y - prev.y);
    X = Vector("X", start, start, INT_MAX);
    X.end.offset(Y.j, -Y.i);
    X.refresh();
    shape = find_shape(segments_j, points, start, size, 'r', 'r', X, Y);
    /* skip NULL shapes */
    if(shape) {
        polygon.shape = *shape;
        polygon.perimeter = find_perimeter(*shape, points);
        tessellations.push_back(polygon);
    }
    /* find the fourth shape using segments[j][1] */
    start = points[segments_j[j][1]];
    prev = points[segments_j[j][0]];
    Y = Vector("Y", prev, start, INT_MAX);
    Y.offset(start.x - prev.x, start.y - prev.y);
    X = Vector("X", start, start, INT_MAX);
    X.end.offset(Y.j, -Y.i);
    X.refresh();
    shape = find_shape(segments_j, points, start, size, 'r', 'r', X, Y);
    /* skip NULL shapes */
    if(shape) {
        polygon.shape = *shape;
        polygon.perimeter = find_perimeter(*shape, points);
        tessellations.push_back(polygon);
    }
    return tessellations;
}

/* finds a shape given the index of the starting point and direction for initializing and adding to the path */
vector<int> *find_shape(vector<int *> segments, Point *points, Point start, int size, char init, char add, Vector X, Vector Y)
{
    Point prev = start;
    Point root = start;
    vector<int *> edges;
    vector<int> *shape = new vector<int>;
    vector<int> path; //contains the path of nodes to add as a polygon
    int j = 0;
    int complete = 0;
    /* find the initial cluster of edges */
    edges = edge_search(segments, start.index, points, size);
    /* find the initial direction */
    path.push_back(edges[0][0]);
    path = init_path(path, edges, points, X, Y, init, 0);
    if(path.size() == 0) {
        printf("ERROR! Path from %d could not be initialized...\n", start.index);
        return NULL;
    }
    j = 1;
    complete = 0;
    while(!complete) {
        prev = start;
        start = points[path[j]];
        /* initialize axis vectors */ //POSSIBLY REVERSE PREV AND START???
        Y = Vector("Y", prev, start, INT_MAX);
        Y.offset(start.x - prev.x, start.y - prev.y);
        X = Vector("X", start, start, INT_MAX);
        X.end.offset(Y.j, -Y.i);
        X.refresh();
        /* find the initial tessellations of edges */
        edges = edge_search(segments, start.index, points, size);
        path = add_path(path, edges, points, X, Y, add);
        if(path.size() == 0) {
            printf("ERROR! Path from %d could not be added...\n", start.index);
            return NULL;
        }
        if(path[j + 1] == path[0]) {
            complete = 1;
        }
        j++;
    }
    /* return the found shape */
    *shape = path;
    return shape;
}

/* initializes the path from a vertex */
vector<int> init_path(vector<int> path, vector<int *> edges, Point *points, Vector X, Vector Y, char type, int branch)
{
    Vector E; //edge vector
    double *X_flags = new double [edges.size()];
    double *Y_flags = new double [edges.size()];
    double dot = 0.0;
    int *quads = new int [edges.size()];
    int i = 0;
    int curr = -1;
    int init = 0;
    int count = -1;
    /* calculate and set flags for each edge */
    for(i = 0; i < edges.size(); i++) {
        /* set edge vector */
        E = Vector("E", points[edges[i][0]], points[edges[i][1]]);
        X_flags[i] = dot_product(E, X);
        Y_flags[i] = dot_product(E, Y);
        //printf("\n");
        //print_v(E);
        //print_v(X);
        //print_v(Y);
        //printf("DOTX: %0.2lf, ", X_flags[i]);
        //printf("DOTY: %0.2lf, ", Y_flags[i]);
        if(dot_product(E, X) > 0.000001) {
            if(dot_product(E, Y) > 0.000001) {
                quads[i] = 1;
            }
            else {
                quads[i] = 4;
            }
        }
        else if(dot_product(E, X) < -0.000001) {
            if(dot_product(E, Y) > 0.000001) {
                quads[i] = 2;
            }
            else {
                quads[i] = 3;
            }
        }
        else {
            if(dot_product(E, Y) > 0.000001) {
                quads[i] = 5;
            }
            else {
                quads[i] = 6;
            }
        }
        //printf("QUAD: %d\n", quads[i]);
    }
    /* check for which direction to initialize to */
    switch(type) {
    case 'r':
        /* first check for an edge in both quadrants 1 and 2 */
        for(i = 0; i < edges.size(); i++) {
            if(quads[i] == 5) {
                count++;
                curr = i;
                break;
            }
        }
        /* check for an edge in quadrant 1 */
        if(count != branch) {
            for(i = 0; i < edges.size(); i++) {
                if(quads[i] == 1) {
                    count++;
                    if(!init) {
                        curr = i;
                        init = 1;
                    }
                    else if(Y_flags[i] > Y_flags[curr]) {
                        curr = i;
                    }
                    if(count == branch) {
                        break;
                    }
                }
            }
        }
        /* check for an edge in quadrant 4 */
        if(count != branch) {
            for(i = 0; i < edges.size(); i++) {
                if(quads[i] == 4) {
                    count++;
                    if(!init) {
                        curr = i;
                        init = 1;
                    }
                    else if(Y_flags[i] > Y_flags[curr]) {
                        curr = i;
                    }
                    if(count == branch) {
                        break;
                    }
                }
            }
        }
        /* check for an edge in both quadrants 3 and 4 */
        if(count != branch) {
            for(i = 0; i < edges.size(); i++) {
                if(quads[i] == 6) {
                    count++;
                    curr = i;
                    break;
                }
            }
        }
        break;
    case 'l':
        /* first check for an edge in both quadrants 1 and 2 */
        for(i = 0; i < edges.size(); i++) {
            if(quads[i] == 5) {
                count++;
                curr = i;
                break;
            }
        }
        /* check for an edge in quadrant 2 */
        if(count != branch) {
            for(i = 0; i < edges.size(); i++) {
                if(quads[i] == 2) {
                    count++;
                    if(!init) {
                        curr = i;
                        init = 1;
                    }
                    else if(Y_flags[i] > Y_flags[curr]) {
                        curr = i;
                    }
                    if(count == branch) {
                        break;
                    }
                }
            }
        }
        /* check for an edge in quadrant 3 */
        if(count != branch) {
            for(i = 0; i < edges.size(); i++) {
                if(quads[i] == 3) {
                    count++;
                    if(!init) {
                        curr = i;
                        init = 1;
                    }
                    else if(Y_flags[i] > Y_flags[curr]) {
                        curr = i;
                    }
                    if(count == branch) {
                        break;
                    }
                }
            }
        }
        /* check for an edge in both quadrants 3 and 4 */
        if(count != branch) {
            for(i = 0; i < edges.size(); i++) {
                if(quads[i] == 6) {
                    count++;
                    curr = i;
                    break;
                }
            }
        }
        break;
    default:
        printf("\nDIRECTION NOT SET! Exiting program...\n\n");
        exit(EXIT_FAILURE);
    }
    /* use curr to create the link in the path */
    if((curr > -1) && (count == branch)) {
        path.push_back(edges[curr][1]);
        if(duplicate_search(path)) {
            //printf("duplicate found: %d\n", edges[curr][1]);
            path.clear();
            return path;
        }
    }
    else {
        path.clear();
    }
    return path;
}

/* adds to the path from a vertex */
vector<int> add_path(vector<int> path, vector<int *> edges, Point *points, Vector X, Vector Y, char type)
{
    Vector E; //edge vector
    double *X_flags = new double [edges.size()];
    double *Y_flags = new double [edges.size()];
    double dot = 0.0;
    int *quads = new int [edges.size()];
    int i = 0;
    int curr = -1;
    int init = 0;
    int found = 0;
    /* initialize the quadrant flags for each edge */
    for(i = 0; i < edges.size(); i++) {
        X_flags[i] = 0.0;
        Y_flags[i] = 0.0;
    }
    /* calculate and set flags for each edge */
    for(i = 0; i < edges.size(); i++) {
        /* set edge vector */ // POSSIBLY SET I AND J TO 1???
        E = Vector("E", points[edges[i][0]], points[edges[i][1]]);
        X_flags[i] = dot_product(E, X);
        Y_flags[i] = dot_product(E, Y);
        if(dot_product(E, X) > 0) {
            if(dot_product(E, Y) > 0) {
                quads[i] = 1;
            }
            else {
                quads[i] = 4;
            }
        }
        else if(dot_product(E, X) < 0) {
            if(dot_product(E, Y) > 0) {
                quads[i] = 2;
            }
            else {
                quads[i] = 3;
            }
        }
        else {
            if(dot_product(E, Y) > 0) {
                quads[i] = 5;
            }
            /* disregaurd segment if it is directly behind */
            else {
                quads[i] = 0;
            }
        }
    }
    /* check for which direction to add */
    switch(type) {
    case 'r':
        /* first check for an edge in quadrant 4 */
        for(i = 0; i < edges.size(); i++) {
            if(quads[i] == 4) {
                found = 1;
                if(!init) {
                    curr = i;
                    init = 1;
                }
                /* find the smallest Y value */
                else if(Y_flags[i] < Y_flags[curr]) {
                    curr = i;
                }
            }
        }
        /* now check for an edge in quadrant 1 */
        if(!found) {
            for(i = 0; i < edges.size(); i++) {
                if(quads[i] == 1) {
                    found = 1;
                    if(!init) {
                        curr = i;
                        init = 1;
                    }
                    /* find the smallest Y value */
                    else if(Y_flags[i] < Y_flags[curr]) {
                        curr = i;
                    }
                }
            }
        }
        /* now check if it is in both quadrant 1 and 2 */
        if(!found) {
            for(i = 0; i < edges.size(); i++) {
                if(quads[i] == 5) {
                    curr = i;
                    break;
                }
            }
        }
        /* now check for an edge in quadrant 2 */
        if(!found) {
            for(i = 0; i < edges.size(); i++) {
                if(quads[i] == 2) {
                    found = 1;
                    if(!init) {
                        curr = i;
                        init = 1;
                    }
                    /* find the greatest Y value */
                    else if(Y_flags[i] > Y_flags[curr]) {
                        curr = i;
                    }
                }
            }
        }
        /* finally check for an edge in quadrant 3 */
        if(!found) {
            for(i = 0; i < edges.size(); i++) {
                if(quads[i] == 3) {
                    found = 1;
                    if(!init) {
                        curr = i;
                        init = 1;
                    }
                    /* find the greatest Y value */
                    else if(Y_flags[i] > Y_flags[curr]) {
                        curr = i;
                    }
                }
            }
        }
        break;
    case 'l':
        /* first check for an edge in quadrant 3 */
        for(i = 0; i < edges.size(); i++) {
            if(quads[i] == 3) {
                found = 1;
                if(!init) {
                    curr = i;
                    init = 1;
                }
                /* find the smallest Y value */
                else if(Y_flags[i] < Y_flags[curr]) {
                    curr = i;
                }
            }
        }
        /* now check for an edge in quadrant 2 */
        if(!found) {
            for(i = 0; i < edges.size(); i++) {
                if(quads[i] == 2) {
                    found = 1;
                    if(!init) {
                        curr = i;
                        init = 1;
                    }
                    /* find the smallest Y value */
                    else if(Y_flags[i] < Y_flags[curr]) {
                        curr = i;
                    }
                }
            }
        }
        /* now check if it is in both quadrant 1 and 2 */
        if(!found) {
            for(i = 0; i < edges.size(); i++) {
                if(quads[i] == 5) {
                    curr = i;
                    break;
                }
            }
        }
        /* now check for an edge in quadrant 1 */
        if(!found) {
            for(i = 0; i < edges.size(); i++) {
                if(quads[i] == 1) {
                    found = 1;
                    if(!init) {
                        curr = i;
                        init = 1;
                    }
                    /* find the greatest Y value */
                    else if(Y_flags[i] > Y_flags[curr]) {
                        curr = i;
                    }
                }
            }
        }
        /* finally check for an edge in quadrant 4 */
        if(!found) {
            for(i = 0; i < edges.size(); i++) {
                if(quads[i] == 4) {
                    found = 1;
                    if(!init) {
                        curr = i;
                        init = 1;
                    }
                    /* find the greatest Y value */
                    else if(Y_flags[i] > Y_flags[curr]) {
                        curr = i;
                    }
                }
            }
        }
        break;
    default:
        printf("\nError: direction not set\n\n");
        exit(EXIT_FAILURE);
    }
    /* use curr to create the link in the path */
    if(curr > -1) {
        path.push_back(edges[curr][1]);
        /* if path is not a viable loop */
        if(duplicate_search(path)) {
            //printf("duplicate found: %d\n", points[edges[curr][1]].index);
            path.clear();
            return path;
        }
    }
    else {
        path.clear();
    }
    return path;
}

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

/* return the index of the matching element from a vector of segments */
int index_match(vector<int *> segments, int vertex)
{
    /* check if the segment is found */
    for(int i = 0; i < segments.size(); i++) {
        if(segments[i][0] == vertex) {
            return i;
        }
        else if(segments[i][1] == vertex) {
            return i;
        }
    }
    return -1;
}

/* searches through a vector of polygons for a matching vertex */
int polygons_search(vector<vector<int> > polygons, int vertex)
{
    int d = 0;
    /* check if the vertex is found */
    for(int i = 0; i < polygons.size(); i++) {
        if(find(polygons[i].begin(), polygons[i].end(), vertex) != polygons[i].end()) {
            return i;
        }
    }
    return -1;
}

/* returns the perimeter for the input shape */
double find_perimeter(vector<int> shape, Point *points)
{
    int i = 0;
    double sum = 0.0;

    /* shape contains indices to points array */
    for(i = 0; i < shape.size() - 1; i++) {
        sum += distance_p(points[shape[i]], points[shape[i + 1]]);
    }
    return sum;
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

/* searches through a vector of shapes for a duplicate */
vector<struct polygon_t> delete_duplicates(vector<struct polygon_t> polygons)
{
    vector<struct polygon_t> tmp = polygons;
    int i = 0;
    int j = 0;
    int k = 0;
    int n = 0;

    /* remove incorrect shapes */
    for(i = 0; i < tmp.size(); i++) {
        if((tmp[i]).shape.size() <= 3) {
            tmp.erase(tmp.begin() + i, tmp.begin() + i + 1);
            polygons.erase(polygons.begin() + i);
            i--;
        }
    }
    /* remove the loop for each shape */
    for(i = 0; i < tmp.size(); i++) {
        (tmp[i]).shape.erase((tmp[i]).shape.begin());
    }
    /* sort each entry */
    for(i = 0; i < tmp.size(); i++) {
        sort((tmp[i]).shape.begin(), (tmp[i]).shape.end());
        //printf("new: ");
        for(j = 0; j < (tmp[i]).shape.size(); j++) {
            //printf("%d ", tmp[i].shape[j]);
        }
        //printf("\n");
    }
    /* check for an equal entry
    */
    for(i = 0; i < tmp.size(); i++) {
        for(j = 0; j < tmp.size(); j++) {
            if(((tmp[i]).shape.size() == (tmp[j]).shape.size()) && (i != j)) {
                for(k = 0; k < (tmp[i]).shape.size(); k++) {
                    if((tmp[i]).shape[k] != (tmp[j]).shape[k]) {
                        break;
                    }
                }
                if(k == (tmp[i]).shape.size()) {
                    //printf("%d. deleting: ", j);
                    for(k = 0; k < (tmp[j]).shape.size(); k++) {
                        //printf("%d ", tmp[j].shape[k]);
                    }
                    //printf("\n");
                    tmp.erase(tmp.begin() + j);
                    polygons.erase(polygons.begin() + j);
                    j--;
                    i--;
                }
            }
        }
    }
    return polygons;
}

/* returns the optimized polygons */
vector<struct polygon_t> optimize_polygons(vector<struct polygon_t> polygons, vector<int *> *segments, Point *points, int size)
{
    Vector V;
    Vector T1;
    Vector T2;
    Point *curr = new Point [size];
    Point best;
    Point start;
    Point end;
    Point prev;
    Point center;
    double sum_x = 0.0;
    double sum_y = 0.0;
    vector<int *> edges;
    int *tmp;
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;
    int n = 0;
    /* initialization of arrays */
    for(i = 0; i < size; i++) {
        curr[i].x = DBL_MAX;
        curr[i].y = DBL_MAX;
        curr[i].tao_distance = DBL_MAX;
        curr[i].index = INT_MAX;
    }
    best.x = DBL_MAX;
    best.y = DBL_MAX;
    best.tao_distance = DBL_MAX;
    best.index = INT_MAX;
    /* make sure each segment branch is optimal */
    for(i = 0; i < polygons.size(); i++) {
        start = Point(points[(polygons[i]).shape[0]]);
        prev = Point(start);
        best = Point(start);
        /* initializing vector T1 */
        T1 = Vector("T1", center, start, INT_MAX);
        T1.offset(start.x - center.x, start.y - center.y);
        /* loops through all segments of the polygon */
        for(j = 0; j < (polygons[i]).shape.size() - 2; j++) {
            k = (polygons[i]).shape[j];
            l = (polygons[i]).shape[j + 1];
            start = Point(points[k]);
            end = Point(points[l]);
            /* initializing vector T1 */
            T1 = Vector("T1", start, end, INT_MAX);
            T1.offset(end.x - start.x, end.y - start.y);
            /* refreshing best index */
            best.tao_distance = DBL_MAX;
            /* initializing vector T2 */
            T2 = Vector("T2", start, start);
            /* initializing vector V */
            T2 = Vector("V", start, start);
            /* loops through all possible indices from start */
            for(k = 0; k < size; k++) {
                /* initializing vector V */
                V.end = Point(points[k]);
                V.refresh();
                /* initializing vector T2 */
                T2 = Vector("T2", T2.start, V.end, INT_MAX);
                T2.offset(V.end.x - T2.start.x, V.end.y - T2.start.y);
                /* initializing tao, theta, and curvature */
                curr[k].tao = (dot_product(T1, T2)); //length of T1 and T2 is always 1
                if(curr[k].tao <= -1.0) {
                    curr[k].tao = -1.0;
                }
                else if(curr[k].tao >= 1.0) {
                    curr[k].tao = 1.0;
                }
                curr[k].x = V.end.x;
                curr[k].y = V.end.y;
                curr[k].index = V.end.index;
                curr[k].theta = angle_t(curr[k].tao);
                curr[k].curvature = calculate_curvature(T1, T2, curr[k].tao);
                curr[k].tao_distance = tao_distance(V, curr[k].curvature, curr[k].theta);
                V.end.tao_distance = curr[k].tao_distance;
            }
            /* find point with the lowest tao-distance */
            for(k = 0; k < size; k++) {
                if(best.tao_distance > curr[k].tao_distance) {
                    best = Point(curr[k]);
                    n = k;
                }
            }
            /* record path */
            for(k = 0; k < (*segments).size(); k++) {
                if(segment_match(*segments, (polygons[i]).shape[j], n) == -1) {
                    tmp = new int [2];
                    tmp[0] = (polygons[i]).shape[j];
                    tmp[1] = n;
                    (*segments).push_back(tmp);
                    //printf("added: <%d,%d>\n", points[(polygons[i]).shape[j]].index, points[n].index);
                    break;
                }
            }
        }
    }
    /* finds the "positive" optimal segments */
    for(i = 0; i < polygons.size(); i++) {
        sum_x = 0;
        sum_y = 0;
        /* calculate average point */
        for(j = 0; j < (polygons[i]).shape.size(); j++) {
            sum_x += points[polygons[i].shape[j]].x;
            sum_y += points[polygons[i].shape[j]].y;
        }
        center.x = sum_x / (polygons[i]).shape.size();
        center.y = sum_y / (polygons[i]).shape.size();
        start = Point(points[(polygons[i]).shape[0]]);
        prev = Point(start);
        best = Point(start);
        /* initializing vector T1 */
        T1 = Vector("T1", center, start, INT_MAX);
        T1.offset(start.x - center.x, start.y - center.y); // POSSIBLY REVERSED???
        /* loops through all indices of a single polygon */
        for(j = 0; j < (polygons[i]).shape.size(); j++) {
            k = (polygons[i]).shape[j];
            start = Point(points[k]);
            /* initializing vector T1 */
            T1 = Vector("T1", center, start, INT_MAX);
            T1.offset(start.x - center.x, start.y - center.y); // POSSIBLY REVERSED???
            /* refreshing best index */
            best.tao_distance = DBL_MAX;
            /* initializing vector T2 */
            T2 = Vector("T2", start, start);
            /* initializing vector V */
            V = Vector("V", start, start);
            /* loops through all possible indices from start */
            for(k = 0; k < size; k++) {
                /* initializing vector V */
                V.end = Point(points[k]);
                V.refresh();
                /* initializing vector T2 */
                T2 = Vector("T2", T2.start, V.end, INT_MAX);
                T2.offset(V.end.x - T2.start.x, V.end.y - T2.start.y);
                /* initializing tao, theta, and curvature */
                curr[k].tao = (dot_product(T1, T2)); //length of T1 and T2 is always 1
                if(curr[k].tao <= -1.0) {
                    curr[k].tao = -1.0;
                }
                else if(curr[k].tao >= 1.0) {
                    curr[k].tao = 1.0;
                }
                curr[k].x = V.end.x;
                curr[k].y = V.end.y;
                curr[k].index = V.end.index;
                curr[k].theta = angle_t(curr[k].tao);
                curr[k].curvature = calculate_curvature(T1, T2, curr[k].tao);
                curr[k].tao_distance = tao_distance(V, curr[k].curvature, curr[k].theta);
                V.end.tao_distance = curr[k].tao_distance;
            }
            /* find point with the lowest tao-distance */
            for(k = 0; k < size; k++) {
                if(best.tao_distance > curr[k].tao_distance) {
                    best = Point(curr[k]);
                    n = k;
                }
            }
            /* record path */
            for(k = 0; k < (*segments).size(); k++) {
                if(segment_match(*segments, (polygons[i]).shape[j], n) == -1) {
                    tmp = new int [2];
                    tmp[0] = (polygons[i]).shape[j];
                    tmp[1] = n;
                    (*segments).push_back(tmp);
                    //printf("added: <%d,%d>\n", points[(polygons[i]).shape[j]].index, points[n].index);
                    break;
                }
            }
        }
    }
    /* finds the "negative" optimal segments */
    for(i = 0; i < polygons.size(); i++) {
        sum_x = 0;
        sum_y = 0;
        /* calculate average point */
        for(j = 0; j < (polygons[i]).shape.size(); j++) {
            sum_x += points[polygons[i].shape[j]].x;
            sum_y += points[polygons[i].shape[j]].y;
        }
        center.x = sum_x / (polygons[i]).shape.size();
        center.y = sum_y / (polygons[i]).shape.size();
        start = Point(points[(polygons[i]).shape[0]]);
        prev = Point(start);
        best = Point(start);
        /* loops through all indices of a single polygon */
        for(j = 0; j < (polygons[i]).shape.size(); j++) {
            k = (polygons[i]).shape[j];
            start = Point(points[k]);
            /* initializing vector T1 */
            T1 = Vector("T1", center, start, INT_MAX);
            T1.offset(-(start.x - center.x), -(start.y - center.y));//POSSIBLY NEED TO SWITCH???
            /* refreshing best index */
            best.tao_distance = DBL_MAX;
            /* initializing vector T2 */
            T2 = Vector("T2", start, start);
            /* initializing vector V */
            V = Vector("V", start, start);
            /* loops through all possible indices from start */
            for(k = 0; k < size; k++) {
                /* initializing vector V */
                V.end = Point(points[k]);
                V.refresh();
                /* initializing vector T2 */
                T2 = Vector("T2", T2.start, V.end, INT_MAX);
                T2.offset(V.end.x - T2.start.x, V.end.y - T2.start.y);
                /* initializing tao, theta, and curvature */
                curr[k].tao = (dot_product(T1, T2)); //length of T1 and T2 is always 1
                if(curr[k].tao <= -1.0) {
                    curr[k].tao = -1.0;
                }
                else if(curr[k].tao >= 1.0) {
                    curr[k].tao = 1.0;
                }
                curr[k].x = V.end.x;
                curr[k].y = V.end.y;
                curr[k].index = V.end.index;
                curr[k].theta = angle_t(curr[k].tao);
                curr[k].curvature = calculate_curvature(T1, T2, curr[k].tao);
                curr[k].tao_distance = tao_distance(V, curr[k].curvature, curr[k].theta);
                V.end.tao_distance = curr[k].tao_distance;
            }
            /* find point with the lowest tao-distance */
            for(k = 0; k < size; k++) {
                if(best.tao_distance > curr[k].tao_distance) {
                    best = Point(curr[k]);
                    n = k;
                }
            }
            /* record path */
            for(k = 0; k < (*segments).size(); k++) {
                if(segment_match(*segments, (polygons[i]).shape[j], n) == -1) {
                    tmp = new int [2];
                    tmp[0] = (polygons[i]).shape[j];
                    tmp[1] = n;
                    (*segments).push_back(tmp);
                    //printf("added: <%d,%d>\n", points[(polygons[i]).shape[j]].index, points[n].index);
                    break;
                }
            }
        }
    }
    /* get rid of single-point segments */
    for(i = 0; i < segments->size(); i++) {
        if((*segments)[i][0] == (*segments)[i][1]) {
            segments->erase(segments->begin() + i);
            i--;
        }
    }
    return polygons;
}

/* returns the polygons with all crosses removed */
void remove_crosses(vector<int *> *segments, Point *points, int size)
{
    struct polygon_t tmp_polygon;
    Point p1;
    Point p2;
    Point p3;
    Point p4;
    Point tmp;
    double y = 0.0;
    double x = 0.0;
    double m1 = 0.0;
    double m2 = 0.0;
    double b1 = 0.0;
    double b2 = 0.0;
    double r = DBL_MAX;
    double l = DBL_MIN;
    double pos_epsilon = 0.000001;
    double neg_epsilon = -0.000001;
    double sum_i = 0.0;
    double sum_j = 0.0;
    double min = DBL_MAX;
    vector<struct polygon_t> cross;
    vector<int *> edges;
    vector<int> tmp_shape;
    int erase = 0;
    int i = 0;
    int j = 0;
    int k = 0;
    int n = 0;
    int m = 0;
    /* loops through all segments */
    for(i = 0; i < segments->size(); i++) {
        /* loops through all segments */
        for(j = 0; j < segments->size(); j++) {
            if(i == j) {
                continue;
            }
            p1 = points[(*segments)[i][0]];
            p2 = points[(*segments)[i][1]];
            p3 = points[(*segments)[j][0]];
            p4 = points[(*segments)[j][1]];
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
                continue;
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
                continue;
            }
            /* p3---------p4  p1---------p2 */
            if(l >= r) {
                continue;
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
            /* erase the segment if it is an overlap */
            if((b2 == b1) && (m1 == m2)) {
                if(p3.x >= p2.x) {
                    break;
                }
                else if(p1.x >= p4.x) {
                    break;
                }
            }
            if(x <= r && x >= l) {
                /* store shapes of cross */
                printf("CROSS (%d,%d):(%d,%d)\n", points[(*segments)[i][0]].index, points[(*segments)[i][1]].index, points[(*segments)[j][0]].index, points[(*segments)[j][1]].index);
                cross = tessellate_cross(*segments, i, j, points, size);
                for(k = 0; k < cross.size(); k++) {
                    printf("%d: ", k);
                    for(n = 0; n < cross[k].shape.size() - 1; n++) {
                        printf("%d->", cross[k].shape[n]);
                    }
                    printf("%d = %lf\n", cross[k].shape[n], cross[k].perimeter);
                }
                /* find the first smallest shape */
                min = DBL_MAX;
                for(k = 0; k < cross.size(); k++) {
                    if(min >= cross[k].perimeter) {
                        n = k;
                        min = cross[k].perimeter;
                    }
                }
                min = DBL_MAX;
                /* find the second smallest shape */
                for(k = 0; k < cross.size(); k++) {
                    if((min >= cross[k].perimeter) && (k != n)) {
                        m = k;
                        min = cross[k].perimeter;
                    }
                }
                printf("CHOSEN SHAPES: %d, %d\n", n, m);
                printf("Comparing shape perimeters...");
                /* erase the segment with the least shape perimeter */
                if((cross[n].perimeter - cross[m].perimeter) > pos_epsilon) { //erase polygon m
                    printf("erasing %d, ", m);
                    if(edge_match(cross[m], (*segments)[j]) > -1) {
                        printf("(%d, %d)!\n", points[(*segments)[j][0]].index, points[(*segments)[j][1]].index);
                        if(i > j) {
                            i--;
                        }
                        segments->erase(segments->begin() + j);
                        j--;
                    }
                    else {
                        printf("(%d, %d)!\n", points[(*segments)[i][0]].index, points[(*segments)[i][1]].index);
                        if(i < j) {
                            j--;
                        }
                        segments->erase(segments->begin() + i);
                        i--;
                    }
                }
                else if((cross[n].perimeter - cross[m].perimeter) < neg_epsilon) { //erase polygon n
                    printf("erasing %d, ", n);
                    if(edge_match(cross[n], (*segments)[j]) > -1) {
                        printf("(%d, %d)!\n", points[(*segments)[j][0]].index, points[(*segments)[j][1]].index);
                        if(i > j) {
                            i--;
                        }
                        segments->erase(segments->begin() + j);
                        j--;
                    }
                    else {
                        printf("(%d, %d)!\n", points[(*segments)[i][0]].index, points[(*segments)[i][1]].index);
                        if(i < j) {
                            j--;
                        }
                        segments->erase(segments->begin() + i);
                        i--;
                    }
                }
                else { //continue to next phase
                    printf("equal!\n");
                    printf("Comparing segment lengths...");
                    /* compare the segment lengths */
                    if((distance_p(p1, p2) - distance_p(p3, p4)) < neg_epsilon) { //erase the j-segment, i.e. j.val > i.val
                        printf("erasing j-segment!\n");
                        if(i > j) {
                            i--;
                        }
                        segments->erase(segments->begin() + j);
                        j--;
                    }
                    else if((distance_p(p1, p2) - distance_p(p3, p4)) > pos_epsilon) { //erase the i-segment, i.e. i.val > j.val
                        printf("erasing i-segment!\n");
                        if(i < j) {
                            j--;
                        }
                        segments->erase(segments->begin() + i);
                        i--;
                    }
                    else { //continue to next phase
                        printf("equal!\n");
                        printf("Comparing perimeter sums...");
                        sum_i = 0.0;
                        sum_j = 0.0;
                        /* calcualte the perimeter sum for i-segment */
                        edges = edge_search(*segments, points[(*segments)[i][0]].index, points, size);
                        for(k = 0; k < edges.size(); k++) {
                            sum_i += distance_p(points[(*segments)[i][0]], points[edges[k][1]]);
                        }
                        edges = edge_search(*segments, points[(*segments)[i][1]].index, points, size);
                        for(k = 0; k < edges.size(); k++) {
                            sum_i += distance_p(points[(*segments)[i][1]], points[edges[k][1]]);
                        }
                        /* calcualte the perimeter sum for j-segment */
                        edges = edge_search(*segments, points[(*segments)[j][0]].index, points, size);
                        for(k = 0; k < edges.size(); k++) {
                            sum_j += distance_p(points[(*segments)[j][0]], points[edges[k][1]]);
                        }
                        edges = edge_search(*segments, points[(*segments)[j][1]].index, points, size);
                        for(k = 0; k < edges.size(); k++) {
                            sum_j += distance_p(points[(*segments)[j][1]], points[edges[k][1]]);
                        }
                        if((sum_i - sum_j) < neg_epsilon) { //erase the j-segment, i.e. i.val < j.val
                            printf("erasing j-segment!\n");
                            if(i > j) {
                                i--;
                            }
                            segments->erase(segments->begin() + j);
                            j--;
                        }
                        else if((sum_i - sum_j) > pos_epsilon) { //erase the i-segment, i.e. j.val < i.val
                            printf("erasing i-segment!\n");
                            if(i < j) {
                                j--;
                            }
                            segments->erase(segments->begin() + i);
                            i--;
                        }
                        else { //erase both segments
                            printf("erasing i-segment and j-segment!\n");
                            if(i < j) {
                                segments->erase(segments->begin() + j);
                                j--;
                                segments->erase(segments->begin() + i);
                            }
                            else {
                                segments->erase(segments->begin() + i);
                                i--;
                                segments->erase(segments->begin() + j);
                            }
                            i--;
                            j--;
                        }
                    }
                }
            }
        }
    }
}

void finalize_segments(vector<int *> *segments, Point *points, int size)
{
    vector<int *> tmp_segments;
    vector<int *> added_segments;
    int *pushed_segment = NULL;
    vector<int> remove;
    int push = 1;
    int keep = 1;
    int i = 0;
    int j = 0;
    int k = 0;
    int n = 0;
    printf("\nFINALIZING SEGMENTS!!!\n\n");
    /* adding border segments */
    for(i = 0; i < size; i++) {
        /* goes through all prospective points */
        for(j = 0; j < size; j++) {
            push = 1;
            if(i == j) {
                continue;
            }
            /* skip segments that are already recorded */
            if((segment_match(*segments, i, j) > -1) || (segment_match(tmp_segments, i, j) > -1)) {
                continue;
            }
            /* goes through all current segments */
            for(k = 0; k < segments->size(); k++) {
                if(intersection(Vector("V1", points[i], points[j]), Vector("V2", points[(*segments)[k][0]], points[(*segments)[k][1]]))) {
                    push = 0;
                    break;
                }
            }
            /* goes through all temporary segments */
            for(k = 0; k < tmp_segments.size(); k++) {
                if(intersection(Vector("V1", points[i], points[j]), Vector("V2", points[tmp_segments[k][0]], points[tmp_segments[k][1]]))) {
                    if(find(remove.begin(), remove.end(), k) == remove.end()) {
                        remove.push_back(k);
                        push = 0;
                    }
                }
            }
            /* pushes segment to temporary segments */
            if(push) {
                pushed_segment = new int [2];
                pushed_segment[0] = i;
                pushed_segment[1] = j;
                tmp_segments.push_back(pushed_segment);
            }
        }
    }
    /* removes all intersecting temporary segments */
    sort(remove.begin(), remove.end());
    for(k = remove.size() - 1; k > -1; k--) {
        tmp_segments.erase(tmp_segments.begin() + remove[k]);
    }
    /* pushing border segments */
    for(k = 0; k < tmp_segments.size(); k++) {
        segments->push_back(tmp_segments[k]);
    }
    tmp_segments = *segments;
    /* adding edges via removing crosses */
    for(i = 0; i < size; i++) {
        /* goes through all prospective points */
        for(j = 0; j < size; j++) {
            push = 1;
            if(i == j) {
                continue;
            }
            /* skip segments that are already recorded */
            if(segment_match(*segments, i, j) > -1) {
                continue;
            }
            /* goes through all current segments */
            for(k = 0; k < segments->size(); k++) {
                if(intersection(Vector("V1", points[i], points[j]), Vector("V2", points[(*segments)[k][0]], points[(*segments)[k][1]]))) {
                    push = 0;
                    break;
                }
            }
            /* adds segment if it does not intersect current segments */
            if(push) {
                pushed_segment = new int [2];
                pushed_segment[0] = i;
                pushed_segment[1] = j;
                added_segments.push_back(pushed_segment);
            }
        }
    }
    /* return if there are no segments to add */
    if(added_segments.size() == 0) {
        return;
    }
    /* adds one segment at a time and optimizes with remove_crosses */
    tmp_segments = *segments;
    /* push the segment if there is only one */
    if(added_segments.size() == 1) {
        segments->push_back(added_segments[0]);
    }
    else {
        for(i = 0; i < added_segments.size(); i++) {
            keep = 1;
            for(j = 0; j < added_segments.size(); j++) {
                if((i == j) || (segment_match(tmp_segments, added_segments[j][0], added_segments[j][1]) > -1)) {
                    continue;
                }
                tmp_segments.push_back(added_segments[j]);
                remove_crosses(&tmp_segments, points, size);
                /* don't keep the i-segment if it was removed */
                if(segment_match(tmp_segments, added_segments[i][0], added_segments[i][1]) == -1) {
                    keep = 0;
                    break;
                }
                /* remove the j-segment if it wasn't removed */
                if((n = segment_match(tmp_segments, added_segments[j][0], added_segments[j][1])) > -1) {
                    tmp_segments.erase(tmp_segments.begin() + n);
                }
            }
            if(keep) {
                printf("Adding...(%d,%d)\n", points[added_segments[i][0]].index, points[added_segments[i][1]].index);
                segments->push_back(added_segments[i]);
            }
        }
    }
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

struct polygon_t find_shortest_path(vector<struct polygon_t> polygons, Point *points, int size)
{
    vector<struct polygon_t> original = polygons;
    vector<struct polygon_t> base;
    vector<struct polygon_t> bridge;
    struct polygon_t tmp_polygon;
    struct polygon_t shortest_path;
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
    shortest_path = base[1];
    return shortest_path;
}

/* checks which polygon the segments are contained in */
int accept_polygon(struct polygon_t polygon, vector<int *> segments, Point *points) {
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
int smallest_neighbour(vector<struct polygon_t> polygons, struct polygon_t source, int n)
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
vector<int *> disjoint_edges(struct polygon_t A, struct polygon_t B)
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

vector<int *> shared_edges(struct polygon_t A, struct polygon_t B)
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

vector<int> shared_points(struct polygon_t A, struct polygon_t B)
{
    vector<int> shared;
    for(int i = 0; i < A.shape.size() - 1; i++) {
        if(shape_search(B.shape, A.shape[i]) > -1) {
            shared.push_back(i);
        }
    }
    return shared;
}

void visit_polygon(int *visited, struct polygon_t polygon, Point *points)
{
    for(int i = 0; i < polygon.shape.size(); i++) {
        visited[points[polygon.shape[i]].index] = 1;
    }
}

struct polygon_t add_polygons(struct polygon_t A, struct polygon_t B, Point *points)
{
    struct polygon_t C; //A + B
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
struct polygon_t sub_polygons(struct polygon_t A, struct polygon_t B, Point *points)
{
    struct polygon_t C;
    struct polygon_t tmp_polygon;
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
int edge_match(struct polygon_t polygon, int *edge)
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
vector<int *> midpoint_construction(Point *points, int size, FILE *gnu_files[NUM_FILES]) {

    /////////////////////////////////
    // INITIALIZATION OF VARIABLES //
    /////////////////////////////////

    vector<int *> segments;
    int *tmp_segment = new int [2];
    int i = 0;
    int j = 0;

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

    interval /= 16;
    //printf("\nINTERVAL: %lf\n\n", interval);

    /* test all line segments */
    for(i = 0; i < size; i++) {
        for(j = 0; j < size; j++) {
            if(i == j) {
                continue;
            }
            Vector L = Vector("L", points[i], points[j]);
            if(segment_match(segments, i, j) == -1) {
                /* make sure the segment doesn't overlap any other segments */
                bool found_overlap = false;
                for(i = 0; i < segments.size(); i++) {
                    Vector V = Vector("V", points[segments[i][0]], points[segments[i][1]]);
                    if(overlap(L, V)) {
                        found_overlap = true;
                        break;
                    }
                }
                if(found_overlap) {
                    continue;
                }
                /* if the line is valid, add it as a segment */
                if(test_w_segment(L, interval, points, size)) {
                    /* record segment */
                    tmp_segment = new int [2];
                    tmp_segment[0] = i;
                    tmp_segment[1] = j;
                    segments.push_back(tmp_segment);
                    //printf("ADDING: (%d, %d)\n", points[i].index, points[j].index);
                    //segments = remove_crossing_segments(segments, segments.size() - 1, points);
                }
            }
        }
    }

    return segments;
}

/* tests if a weslean segment is valid
 * @param L, the weslean line being tested
 * @param interval, the weslean point generation interval
 * @param points, the array of datapoints
 * @param n, the size of the array
 * @return ..., the line segement to add, either a line or NULL*/
bool test_w_segment(Vector L, double interval, Point *points, int n) {
    int i = 0;
    int j = 0;
    int k = 0;

    //printf("TESTING: (%d, %d)\n", L.start.index, L.end.index);

    /* create initial point set */
    int l = 0;
    int m = 0;
    Point *t_points = new Point [n + 1];
    for(l = 0; l < n; l++) {
        t_points[m++] = points[l];
    }

    /* create w_points */
    vector<Point> *w_points = new vector<Point>;
    *w_points = generate_w_points(*w_points, L, interval);

    /* test all w_points, q = {q_1, q_2, ..., q*} */
    for(j = 0; j < w_points->size(); j++) {
        Point q = Point((*w_points)[j]);
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
                if(!(determinant(S, R) <= 0.000001 && determinant(S, R) >= -0.000001) && ((determinant(Q, R) < 0 && determinant(S, R) < 0) || (determinant(Q, R) > 0 && determinant(S, R) > 0))) {
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
                return false;
            }
        }
    }
    return true;
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
 * @return the updated segments vector
 */
vector<int *> remove_crossing_segments(vector<int *> segments, int s, Point *points) {
    vector<int> crosses;
    Vector L = Vector("L", points[segments[s][0]], points[segments[s][1]]);
    printf("L: (%d, %d)\n", L.start.index, L.end.index);
    int i;

    /* initialize intersections vector */
    for(i = segments.size() - 1; i >= 0; i--) {
        if(i == s) {
            crosses.push_back(i);
            continue;
        }
        Vector V = Vector("V", points[segments[i][0]], points[segments[i][1]]);
        if(intersection(V, L)) {
            crosses.push_back(i);
        }
    }

    /* short circuit */
    if(crosses.size() <= 1) {
        return segments;
    }

    /* find the minimum distance */
    double min = distance_p(points[segments[crosses[0]][0]], points[segments[crosses[0]][1]]);
    for(i = 1; i < crosses.size(); i++) {
        if(min > distance_p(points[segments[crosses[i]][0]], points[segments[crosses[i]][1]])) {
            min = distance_p(points[segments[crosses[i]][0]], points[segments[crosses[i]][1]]);
        }
    }

    int equal_count = 0;
    for(i = 0; i < crosses.size(); i++) {
        double distance = distance_p(points[segments[crosses[i]][0]], points[segments[crosses[i]][1]]);
        if(min == distance) {
            equal_count++;
        }
    }

    /* delete all segments with length less than min */
    for(i = 0; i < crosses.size(); i++) {
        double distance = distance_p(points[segments[crosses[i]][0]], points[segments[crosses[i]][1]]);
        if(min == distance) {
            /* if more than one remaining segment has length min, delete them all */
            if(equal_count > 1) {
                printf("REMOVING: (%d, %d)\n", points[segments[crosses[i]][0]].index, points[segments[crosses[i]][1]].index);
                segments.erase(segments.begin() + crosses[i]);
            }
        }
        if(min < distance) {
            printf("REMOVING: (%d, %d)\n", points[segments[crosses[i]][0]].index, points[segments[crosses[i]][1]].index);
            segments.erase(segments.begin() + crosses[i]);
        }
    }

    return segments;
}

/* returns true if the points overlap, false otherwise */
bool overlap(Vector V1, Vector V2) {
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
    /* find overlap area */
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
    return false;
}
