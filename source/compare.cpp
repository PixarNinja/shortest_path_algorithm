#include "compare.h"

#ifndef COMPARE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string>
#include <vector>

#include "point.h"

using namespace std;

#endif

/* calculates distance given two points
 * @param P1, the start point
 * @param P2, the end point
 * @return the distance between the points
 */
double distance(Point P1, Point P2) {
    return sqrt(pow(P2.x - P1.x, 2) + pow(P2.y - P1.y, 2));
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

int main(int argc, char *argv[])
{
    if(argc != 4) {
        printf("Usage: ./compare [tessellations path] [shortest path] [datapoint path]\n");
        exit(EXIT_FAILURE);
    }
    FILE *tess = fopen(argv[1], "r");
    FILE *path = fopen(argv[2], "r");
    FILE *data = fopen(argv[3], "r");
    if(tess == NULL || path == NULL || data == NULL) {
        printf("File(s) not found. Exiting Program. Good Day.\n");
        exit(EXIT_FAILURE);
    }
    char buf[1024];
    int i = 0;
    int j = 0;
    int p1 = -1;
    int p2 = -1;
    int p3 = -1;
    int p4 = -1;
    int found = 0;

    int path_size = 0;
    while(fgets(buf, 1024, path)) {
        path_size++;
    }
    int tess_size = 0;
    while(fgets(buf, 1024, tess)) {
        tess_size++;
    }

    Point *points;
    int size = 0;
    while(fgets(buf, 1024, data)) {
        size++;
    }
    fclose(data);
    data = fopen(argv[3], "r");
    points = new Point [size];
    while(fscanf(data, "%d: (%lf, %lf)", &points[i].index, &points[i].x, &points[i].y) > 0) {
        i++;
    }
    fclose(data);

    printf("\n... comparing %s and %s\n", argv[1], argv[2]);

    fclose(path);
    path = fopen(argv[2], "r");
    for(i = 0; i < path_size; i++) {
        fscanf(path, "%d %d", &p1, &p2);
        fclose(tess);
        tess = fopen(argv[1], "r");
        found = 0;
        for(j = 0; j < tess_size; j++) {
            fscanf(tess, "%d %d", &p3, &p4);
            if(((p1 == p3) && (p2 == p4)) || ((p2 == p3) && (p1 == p4))) {
                found = 1;
                break;
            }
        }
        if(!found) {
            printf("... PATH NOT FOUND\n\n");
            /* find distances and calculate error */
            fclose(path);
            path = fopen(argv[2], "r");
            double distance_calc = 0.0;
            for(i = 0; i < path_size; i++) {
                fscanf(path, "%d %d", &p1, &p2);
                distance_calc += distance(points[p1], points[p2]);
            }
            break;
        }
    }
    if(found) {
        printf("... PATH FOUND\n\n");
    }

    fclose(path);
    fclose(tess);
    return 0;
}
