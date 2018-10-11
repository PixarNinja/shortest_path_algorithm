/*
 * Brute force function definitions
 *
 * Shortest Path Algorithm
 * Mark Wesley Harris
 * 2018
 */

#include "brute_force.h"

#ifndef BRUTEFORCE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <vector.h>
#include <algorithm>

#include "point.h"
#include "vector.h"

#endif

int *global_shortest;
int global_count = 0;

int main(int argc, char *argv[])
{
    if(argc != 3) {
        printf("\nUsage: ./brute_force [flag] [datapoint path]\n\n");
        exit(EXIT_FAILURE);
    }
    int flag = atoi(argv[1]);
    FILE *data = fopen(argv[2], "r");
    if(!data) {
        printf("\nERROR: unable to open data file. Exiting Program. Good day.\n\n");
        exit(EXIT_FAILURE);
    }
    struct data_t *point;
    struct data_t *shortest;
    char buf[1024];
    int start = 0;
    int n;
    int size = 0;
    int i = 0;

    const char *gnu_path = "./gnu_files/";
    char *name = new char [1024];

    char *tmp = new char [strlen(argv[2]) + 1];
    for(i = 0; i < strlen(argv[2]) + 1; i++)
        tmp[i] = argv[2][i];
    tmp = strtok(tmp, "/");
    while (tmp != NULL) {
        for(i = 0; i < strlen(tmp) + 1; i++)
            name[i] = tmp[i];
        tmp = strtok (NULL, "/");
    }
    name = strtok(name, ".");
    
    while(fgets(buf, 1024, data)) {
        size++;
    }
    n = size - 1;

    global_shortest = new int [size + 1];
    for(i = 0; i < size; i++) {
        global_shortest[i] = size;
    }

    i = 0;
    fclose(data);
    data = fopen(argv[2], "r");
    if(!data) {
        printf("\nERROR: unable to open data file. Exiting Program. Good day.\n\n");
        exit(EXIT_FAILURE);
    }
    point = new data_t [size];
    while(fscanf(data, "%d: (%lf, %lf)", &point[i].p.index, &point[i].p.x, &point[i].p.y) > 0) {
        point[i].curr = DBL_MAX;
        point[i].i = i;
        point[i].v = Vector("V", point[i].p, point[i].p);
        i++;
    }
    fclose(data);
    shortest = point;

    /* switch between distance or curvature */
    char *path;
    FILE *plot;
    switch(flag) {
    case 1:
        path = new char [(strlen(gnu_path) + strlen(name) + strlen("_shortest_dist.gpf") + 1)];
        sprintf(path, "%s%s_shortest_dist.gpf", gnu_path, name);
        plot = fopen(path, "w+");
        if(!plot) {
            printf("\nERROR: unable to open plot file(s). Exiting Program. Good day.\n\n");
            exit(EXIT_FAILURE);
        }

        sprintf(path, "%s%s_shortest_dist.gpf", gnu_path, name);
        /* fills array with permutations */
        permute_dist(point, shortest, start, n);

        printf("\nSHORTEST DISTANCE:\n");
        break;
    case 2:
        path = new char [(strlen(gnu_path) + strlen(name) + strlen("_shortest_curv.gpf") + 1)];
        sprintf(path, "%s%s_shortest_curv.gpf", gnu_path, name);
        plot = fopen(path, "w+");
        if(!plot) {
            printf("\nERROR: unable to open plot file(s). Exiting Program. Good day.\n\n");
            exit(EXIT_FAILURE);
        }
        /* fills array with permutations */
        permute_curv(point, shortest, start, n);

        printf("\nSHORTEST CURVATURE:\n");
        break;
    case 3:
        path = new char [(strlen(gnu_path) + strlen(name) + strlen("_shortest_tao.gpf") + 1)];
        sprintf(path, "%s%s_shortest_tao.gpf", gnu_path, name);
        plot = fopen(path, "w+");
        if(!plot) {
            printf("\nERROR: unable to open plot file(s). Exiting Program. Good day.\n\n");
            exit(EXIT_FAILURE);
        }
        /* fills array with permutations */
        permute_tao(point, shortest, start, n);

        printf("\nSHORTEST TAO DISTANCE:\n");
        break;
    default:
        printf("ERROR: invalid flag, enter 1 for distance, 2 for curvature, and 3 for tao distance\n\n");
        return EXIT_FAILURE;
    }

    path = new char [strlen(argv[2]) + strlen(".shortest_path")];
    char *file_path = argv[2];
    file_path[strlen(argv[2]) - 4] = '\0';
    sprintf(path, "%s.shortest_path", file_path);
    //printf("%s\n", path);
    FILE *output = fopen(path, "w+");

    /* output data */
    printf("\n%d->", global_shortest[0]);
    fprintf(plot, "%lf %lf\n", point[global_shortest[0]].p.x, point[global_shortest[0]].p.y);
    fprintf(output, "%d %d\n", point[global_shortest[0]].p.index, point[global_shortest[1]].p.index);
    for(i = 1; i < n; i++) {
        printf("%d->", point[global_shortest[i]].p.index);
        fprintf(plot, "%lf %lf\n", point[global_shortest[i]].p.x, point[global_shortest[i]].p.y);
        fprintf(output, "%d %d\n", point[global_shortest[i]].p.index, point[global_shortest[i + 1]].p.index);
    }
    printf("%d->%d", point[global_shortest[i]].p.index, point[global_shortest[0]].p.index);
    fprintf(plot, "%lf %lf\n", point[global_shortest[i]].p.x, point[global_shortest[i]].p.y);
    fprintf(plot, "%lf %lf\n", point[global_shortest[0]].p.x, point[global_shortest[0]].p.y);
    fprintf(output, "%d %d\n", point[global_shortest[i]].p.index, point[global_shortest[0]].p.index);
    printf("Total Permutations: %d\n", global_count);
    printf("Value: %lf\n", shortest[0].curr);
    printf("\n");
    fclose(output);

    return 0;
}

/* finds the angle between another vector
 * @param V1, the first vector
 * @param V2, the second vector
 * @return the angle in radians
 */
double angle(Vector V1, Vector V2) {
    return (acos(dot_product(V1, V2) / (V1.length * V2.length)));
}

/* finds the dot product of two vectors
 * @param V1, the first vector
 * @param V2, the second vector
 * @return the dot product of the vectors
 */
double dot_product(Vector V1, Vector V2) {
    return ((V1.i * V2.i) + (V1.j * V2.j));
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

/* returns the determinant between two vectors
 * @param V1, the first vector
 * @param V2, the second vector
 * @return the determinant, i_1*j_2 - j_1*i_2
 */
double determinant(Vector V1, Vector V2) {
    return ((V1.i * V2.j) - (V2.i * V1.j));
}

/* Function to swap values at two pointers */
void swap(struct data_t *x, struct data_t *y)
{
    struct data_t tmp;
    tmp.p.x = x->p.x;
    tmp.p.y = x->p.y;
    tmp.p.index = x->p.index;
    tmp.i = x->i;
    tmp.curr = x->curr;
    x->p.x = y->p.x;
    x->p.y = y->p.y;
    x->p.index = y->p.index;
    x->i = y->i;
    x->curr = y->curr;
    y->p.x = tmp.p.x;
    y->p.y = tmp.p.y;
    y->p.index = tmp.p.index;
    y->i = tmp.i;
    y->curr = tmp.curr;
}

/* calculates all permutations */
void permute_dist(struct data_t *point, struct data_t *shortest, int index, int n)
{
    int i = 0;
    double total = 0;
    double curr = point[0].curr;
    double segment = 0;
    /* base case */
    if(index == n) {
        global_count++;
        /* calculating distance of segments from start to end */
        for(i = 0; i < n; i++) {
            segment = distance(point, i, i + 1);
            total += segment;
        }
        /* calculating the final segment (start and end nodes) */
        segment = distance(point, n, 0);
        total += segment;
        /* checking to see if the most recent total is less than
           the current shortest length */
        if(total < curr) {
            point[0].curr = total;
            global_shortest[0] = shortest[n].i;
            for(i = 0; i <= n; i++) {
                global_shortest[i + 1] = shortest[i].i;
            }
            shortest = point;
            curr = total;
        }
    }
    else {
        for(i = index; i <= n; i++) {
            swap(point + index, point + i);
            permute_dist(point, shortest, index + 1, n);
            swap(point + index, point + i); //backtrack
        }
    }
}

/* calculates all permutations */
void permute_tao(struct data_t *point, struct data_t *shortest, int index, int n)
{
    int i = 0;
    int j = 0;
    double total = 0;
    double curr = point[0].curr;
    double segment = 0;
    vector<int *> path;
    int *tmp;
    /* base case */
    if(index == n) {
        global_count++;
        for(i = 0; i <= n; i++) {
            point[i].v = Vector("V", point[i].p, point[i].p);
        }
        /* calculating tao of segments from start to end */
        for(i = 0; i < n - 1; i++) {
            segment = tao(point, i, i + 1, i + 2, path);
            if(segment < 0) {
                /* segment crosses another in the path */
                return;
            }
            //printf("(%d, %d, %d) = %lf\n", point[i].p.index, point[i + 1].p.index, point[i + 2].p.index, segment);
            total += segment;
            tmp = new int [2];
            tmp[0] = i;
            tmp[1] = i + 1;
            path.push_back(tmp);
        }
        /* calculating the final 2 segments  */
        segment = tao(point, n - 1, n, 0, path);
        if(segment < 0) {
            /* segment crosses another in the path */
            return;
        }
        total += segment;
        tmp = new int [2];
        tmp[0] = n - 1;
        tmp[1] = n;
        path.push_back(tmp);
        //printf("(%d, %d, %d) = %lf\n\n", point[n - 1].p.index, point[n].p.index, point[0].p.index, segment);
        segment = tao(point, n, 0, 1, path);
        if(segment < 0) {
            /* segment crosses another in the path */
            return;
        }
        total += segment;
        //printf("(%d, %d, %d) = %lf\n\n", point[n].p.index, point[0].p.index, point[1].p.index, segment);
        /* checking to see if the most recent total is less than
           the current shortest length */
        if(total < curr) {
            point[0].curr = total;
            global_shortest[0] = shortest[n].i;
            for(i = 0; i <= n; i++) {
                global_shortest[i + 1] = shortest[i].i;
            }
            shortest = point;
            curr = total;
        }
    }
    else {
        for(i = index; i <= n; i++) {
            swap(point + index, point + i);
            permute_tao(point, shortest, index + 1, n);
            swap(point + index, point + i); //backtrack
        }
    }
}

/* calculates all permutations */
void permute_curv(struct data_t *point, struct data_t *shortest, int index, int n)
{
    int i = 0;
    int j = 0;
    double total = 0;
    double curr = point[0].curr;
    double segment = 0;
    vector<int *> path;
    int *tmp;
    /* base case */
    if(index == n) {
        global_count++;
        for(i = 0; i <= n; i++) {
            point[i].v = Vector("V", point[i].p, point[i].p);
        }
        /* calculating curvature of segments from start to end */
        for(i = 0; i < n - 1; i++) {
            segment = curvature(point, i, i + 1, i + 2, path);
            if(segment < 0) {
                /* segment crosses another in the path */
                return;
            }
            //printf("(%d, %d, %d) = %lf\n", point[i].p.index, point[i + 1].p.index, point[i + 2].p.index, segment);
            total += segment;
            tmp = new int [2];
            tmp[0] = i;
            tmp[1] = i + 1;
            path.push_back(tmp);
        }
        /* calculating the final 2 segments  */
        segment = curvature(point, n - 1, n, 0, path);
        if(segment < 0) {
            /* segment crosses another in the path */
            return;
        }
        total += segment;
        tmp = new int [2];
        tmp[0] = n - 1;
        tmp[1] = n;
        path.push_back(tmp);
        //printf("(%d, %d, %d) = %lf\n\n", point[n - 1].p.index, point[n].p.index, point[0].p.index, segment);
        segment = curvature(point, n, 0, 1, path);
        if(segment < 0) {
            /* segment crosses another in the path */
            return;
        }
        total += segment;
        //printf("(%d, %d, %d) = %lf\n\n", point[n].p.index, point[0].p.index, point[1].p.index, segment);
        /* checking to see if the most recent total is less than
           the current shortest length */
        if(total < curr) {
            point[0].curr = total;
            global_shortest[0] = shortest[n].i;
            for(i = 0; i <= n; i++) {
                global_shortest[i + 1] = shortest[i].i;
            }
            shortest = point;
            curr = total;
        }
    }
    else {
        for(i = index; i <= n; i++) {
            swap(point + index, point + i);
            permute_curv(point, shortest, index + 1, n);
            swap(point + index, point + i); //backtrack
        }
    }
}

/* calculates distance given index and structure */
double distance(struct data_t *point, int first, int last)
{
    double x_1 = point[first].p.x;
    double y_1 = point[first].p.y;
    double x_2 = point[last].p.x;
    double y_2 = point[last].p.y;
    double dist = sqrt((x_2 - x_1)*(x_2 - x_1)+(y_2 - y_1)*(y_2 - y_1));
    return dist;
}

/* calculates curvature given index and structure */
double curvature(struct data_t *point, int prev, int curr, int next, vector<int *> path)
{
    struct data_t *p = &point[prev];
    struct data_t *c = &point[curr];
    struct data_t *n = &point[next];

    if(p->v.i == 0 && p->v.j == 0) {
        /* set c */
        c->v.end = Point(n->p);
        c->v.refresh();
        return 0;
    }
    /* check if the segment crosses another in the path */
    for(int *segment : path) {
        //struct data_t *m = &point[index];
        Vector V = Vector("V", point[segment[0]].p, point[segment[1]].p);
        if(intersection(V, p->v)) {
            return -1;
        }
    }

    /* set c */
    c->v.end = Point(n->p);
    c->v.refresh();

    return (angle(p->v, c->v) * 180 / M_PI);
}

/* calculates curvature given index and structure */
double tao(struct data_t *point, int prev, int curr, int next, vector<int *> path)
{
    struct data_t *p = &point[prev];
    struct data_t *c = &point[curr];
    struct data_t *n = &point[next];

    if(p->v.i == 0 && p->v.j == 0) {
        /* set c */
        c->v.end = Point(n->p);
        c->v.refresh();
        return 0;
    }
    /* check if the segment crosses another in the path */
    for(int *segment : path) {
        //struct data_t *m = &point[index];
        Vector V = Vector("V", point[segment[0]].p, point[segment[1]].p);
        if(intersection(V, p->v)) {
            return -1;
        }
    }

    /* set c */
    c->v.end = Point(n->p);
    c->v.refresh();

    return (distance_p(p->p, c->p) * angle(p->v, c->v) / M_PI);
}

/* calculates distance between two Point objects */
double distance_p(Point P1, Point P2) {
    return sqrt(pow(P2.x - P1.x, 2) + pow(P2.y - P1.y, 2));
}

/* calculates factorial of an integer */
int factorial(int n)
{
    int result = 1;
    int i = 2;
    for(; i <= n; i++) {
        result *= i;
    }
    return result;
}
