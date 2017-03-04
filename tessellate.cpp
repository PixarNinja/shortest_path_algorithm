#include <stdio.h>
#include <stdlib.h>
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

#define NUM_FILES 5

using namespace std;

struct point_t {
    double x;
    double y;
    double theta;
    double curvature;
    double tao;
    double tao_distance;
    int index;
};

struct vector_t {
    char *name;
    struct point_t point[2];
    double length;
    double i;
    double j;
};

struct polygon_t {
    vector<int> shape;
    double perimeter;
};

/* global variables */
int permutations = 1;

void construct_segments(vector<int *> *segments, struct point_t *points, struct point_t begin, int n, int size, FILE *gnu_files[NUM_FILES], int *mapped, int **recorded);
void join_vertex(vector<int *> *segments, struct point_t *points, struct point_t begin, int n, int size, FILE *gnu_files[NUM_FILES]);
void join_segment(vector<int *> *segments, struct point_t *points, struct point_t begin, struct point_t end, int n, int m, int size, FILE *gnu_files[NUM_FILES]);
vector<struct polygon_t> construct_polygons(vector<int *> segments, struct point_t *points, int size, FILE *gnu_files[NUM_FILES]);
vector<vector<int> > tesselate(vector<vector<int> > tesselations, vector<int *> segments, struct point_t *points, int size, char init, char add, FILE *gnu_files[NUM_FILES]);
vector<int> init_path(vector<int> path, vector<int *> edges, struct point_t *points, struct vector_t X, struct vector_t Y, char type);
vector<int> add_path(vector<int> path, vector<int *> edges, struct point_t *points, struct vector_t X, struct vector_t Y, char type);
vector<int *> edge_search(vector<int *> segments, int vertex, struct point_t *points, int size);
int index_match(vector<int *> segments, int vertex);
int shape_search(vector<int> shape, int vertex);
int polygons_search(vector<vector<int> > polygons, int vertex);
double find_perimeter(vector<int> shape, struct point_t *points);
int segment_match(vector<int *> segments, int beginning, int end);
int duplicate_search(vector<int> shape);
vector<struct polygon_t> delete_duplicates(vector<struct polygon_t> polygons);
vector<struct polygon_t> optimize_polygons(vector<struct polygon_t> polygons, vector<int *> *segments, struct point_t *points, int size, FILE *gnu_files[NUM_FILES]);
void remove_crosses(vector<int *> *segments, struct point_t *points, int size, FILE *gnu_files[NUM_FILES]);
int point_match(struct point_t *points, int size, int vertex);
double calculate_curvature(struct vector_t T1, struct vector_t T2, double tao);
double angle_t(double tao);
double tao_distance(struct vector_t V, double curvature, double theta);
double angle_v(struct vector_t V1, struct vector_t V2);
double distance_p(struct point_t start, struct point_t end);
double distance_v(struct vector_t V1, struct vector_t V2);
double length_v(struct vector_t V);
double dot_product(struct vector_t V1, struct vector_t V2);
void print_v(struct vector_t V);
void print(struct vector_t V, struct vector_t T1, struct vector_t T2, double curvature, double theta, double tao, double tao_distance);
void memory_error(void);

int main(int argc, char *argv[])
{
    FILE *data;
    FILE *gnu_files[NUM_FILES];
    vector<struct polygon_t> polygons;
    struct point_t *points;
    struct point_t center;
    char buf[1024];
    double range = 0.0;
    double sum_x = 0.0;
    double sum_y = 0.0;
    vector<int *> *segments = new vector<int *> [1];
    vector<int *> edges;
    int **recorded;
    int *mapped;
    int keep_going = 0;
    int size = 0;
    int i = 0;
    int j = 0;
    int k = 0;
    int count = 0;

    if(argc == 1) {
        printf("\n\nPlease enter the path of the .dat file to read from. Exiting Program. Good Day.\n\n");
    }
    gnu_files[0] = fopen ("./gnu_files/commands.tmp", "w+");
    gnu_files[1] = fopen("./gnu_files/datapoints.tmp", "w+");
    gnu_files[2] = fopen("./gnu_files/lines.tmp", "w+");
    gnu_files[3] = fopen("./gnu_files/extrapoints.tmp", "w+");
    gnu_files[4] = fopen("./gnu_files/polygons.tmp", "w+");
    printf("%s\n", argv[1]);
    data = fopen(argv[1], "r");
    while(fgets(buf, 1024, data)) {
        size++;
    }
    fclose(data);
    points = new struct point_t [size];
    mapped = new int [size];
    recorded = new int * [size];
    /* initializing array */
    for(i = 0; i < size; i++) {
        mapped[i] = 0;
        recorded[i] = new int [size];
        for(j = 0; j < size; j++) {
            recorded[i][j] = 0;
        }
    }
    i = 0;
    data = fopen(argv[1], "r");
    while(fscanf(data, "%d: (%lf, %lf)", &points[i].index, &points[i].x, &points[i].y) > 0) {
        if(fabs(points[i].x) > range) {
            range = fabs(points[i].x);
        }
        if(fabs(points[i].y) > range) {
            range = fabs(points[i].y);
        }
        i++;
    }
    fclose(data);
    /* plot datapoints */
    for(i = 0; i < size; i++) {
        fprintf(gnu_files[1], "%lf %lf %d\n", points[i].x, points[i].y, points[i].index);
    }
    /* plot setup */
    fprintf(gnu_files[0], "set xrange [%lf:%lf]\n", -(range + 1), range + 1);
    fprintf(gnu_files[0], "set yrange [%lf:%lf]\n", -(range + 1), range + 1);
    fprintf(gnu_files[0], "set size ratio 1\n");
    fprintf(gnu_files[0], "set grid\n");
    fprintf(gnu_files[0], "set title \"Shortest Path Algorithm\"\n");
    fprintf(gnu_files[0], "set style line 1 lc rgb \"black\" lw 1\n");
    /* runs tao-distance algorithm on data */
    for(i = 0; i < size; i++) {
        for(j = 0; j < size; j++) {
            if(mapped[j] == 0) {
                keep_going = 1;
            }
        }
        if(keep_going == 1) {
            construct_segments(segments, points, points[i], i, size, gnu_files, mapped, recorded);
            keep_going = 0;
        }
        else {
            i = size;
        }
    }
    /* cleanup unmapped points */
    for(i = 0; i < size; i++) {
        if(mapped[i] == 0) {
            join_vertex(segments, points, points[i], i, size, gnu_files);
        }
    }
    /* cleanup points with only one segment */
    for(i = 0; i < size; i++) {
        edges = edge_search(*segments, points[i].index, points, size);
        if(edges.size() == 1) {
            printf("CLEANUP: <%d,%d>\n", points[edges[0][1]].index, points[i].index);
            join_segment(segments, points, points[edges[0][1]], points[i], edges[0][1], i, size, gnu_files);
        }
    }
    /* get rid of crossing lines */
    remove_crosses(segments, points, size, gnu_files);
    /* optimize tesselations */
    do {
        /* find polygons */
        polygons = construct_polygons(*segments, points, size, gnu_files);
        polygons = delete_duplicates(polygons);
        /* ensure each polygon is optimal */
        count = segments->size();
        polygons = optimize_polygons(polygons, segments, points, size, gnu_files);
        /* get rid of crossing lines */
        remove_crosses(segments, points, size, gnu_files);
    } while(count < segments->size());
    /* find polygons again */
    polygons = construct_polygons(*segments, points, size, gnu_files);
    polygons = delete_duplicates(polygons);
    /* print segment information */
    printf("\nSEGMENTS:\n");
    for(i = 0; i < segments->size(); i++) {
        printf("%d: <%d,%d>\n", i, points[(*segments)[i][0]].index, points[(*segments)[i][1]].index);
    }
    /* plot segment information */
    for(i = 0; i < segments->size(); i++) {
        fprintf(gnu_files[2], "%lf %lf %d\n", points[(*segments)[i][0]].x, points[(*segments)[i][0]].y, points[(*segments)[i][0]].index);
        fprintf(gnu_files[2], "%lf %lf %d\n", points[(*segments)[i][1]].x, points[(*segments)[i][1]].y, points[(*segments)[i][1]].index);
        fprintf(gnu_files[2], "\n");
    }
    /* plot perimeter data */
    for(i = 0; i < polygons.size(); i++) {
        sum_x = 0.0;
        sum_y = 0.0;
        for(j = 1; j < (polygons[i]).shape.size(); j++) {
            sum_x += points[(polygons[i]).shape[j]].x;
            sum_y += points[(polygons[i]).shape[j]].y;
        }
        center.x = sum_x / (double)((polygons[i]).shape.size() - 1);
        center.y = sum_y / (double)((polygons[i]).shape.size() - 1);
        fprintf(gnu_files[3], "%lf %lf %0.2lf\n", center.x, center.y, polygons[i].perimeter);
    }
    /* print polygon information */
    printf("\nPOLYGONS:\n");
    for(i = 0; i < polygons.size(); i++) {
        printf("%d: ", i);
        for(j = 0; j < (polygons[i]).shape.size(); j++) {
            printf("%d ", points[(polygons[i]).shape[j]].index);
        }
        printf("= %0.2lf\n", polygons[i].perimeter);
    }
    /* plot */
    fprintf(gnu_files[0], "plot './gnu_files/lines.tmp' using 1:2 with lines ls 1 title \"shortest path\",");
    fprintf(gnu_files[0], "'./gnu_files/datapoints.tmp' using 1:2 with points pt 7 notitle,");
    fprintf(gnu_files[0], "'' using 1:2:3 with labels point pt 7 offset char -1,-1 notitle,");
    fprintf(gnu_files[0], "'./gnu_files/extrapoints.tmp' using 1:2:3 with labels point pt 3 offset char -1,-1 notitle\n");
    printf("\n");
    printf("Total Permutations: %d\n", permutations);
    printf("\n");
    fclose(gnu_files[0]);
    fclose(gnu_files[1]);
    fclose(gnu_files[2]);
    fclose(gnu_files[3]);
    system("gnuplot -persistent ./gnu_files/commands.tmp");
    fclose(gnu_files[4]);
    return 0;
}

/* runs tao-distance algorithm on dataset and generates optimal segments */
void construct_segments(vector<int *> *segments, struct point_t *points, struct point_t begin, int n, int size, FILE *gnu_files[NUM_FILES], int *mapped, int **recorded)
{
    struct vector_t V;
    struct vector_t T1;
    struct vector_t T2;
    struct point_t *curr = new struct point_t [size];
    struct point_t best;
    struct point_t start;
    struct point_t prev;
    struct point_t center;
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
    best.index = INT_MAX;
    /* calculate average point */
    for(i = 0; i < size; i++) {
        sum_x += points[i].x;
        sum_y += points[i].y;
    }
    center.x = sum_x / size;
    center.y = sum_y / size;
    /* plot center point */
    fprintf(gnu_files[3], "%lf %lf %s\n", center.x, center.y, "C");
    printf("begin = %d\n", begin.index);
    printf("\n");
    start.x = begin.x;
    start.y = begin.y;
    start.index = begin.index;
    prev.x = begin.x;
    prev.y = begin.y;
    prev.index = begin.index;
    best.x = begin.x;
    best.y = begin.y;
    best.index = begin.index;
    loop[m++] = n;
    /* initializing vector T1 */
    T1.point[0].x = start.x;
    T1.point[0].y = start.y;
    T1.point[0].index = start.index;
    T1.i = (center.x - start.x) / distance_p(center, start);
    T1.j = (center.y - start.y) / distance_p(center, start);
    T1.point[1].x = start.x + T1.i;
    T1.point[1].y = start.y + T1.j;
    T1.point[1].index = INT_MAX;
    T1.length = length_v(T1);
    /* outer loop, calculates total distance */
    while(visited_count <= total_size) {
        /* store start index in visited-array */
        visited[start.index] = 1;
        i = 0;
        /* refreshing best index */
        best.tao_distance = DBL_MAX;
        best.index = start.index;
        /* initializing vector T2 */
        T2.point[0].x = start.x;
        T2.point[0].y = start.y;
        T2.point[0].index = start.index;
        T2.i = 0;
        T2.j = 0;
        T2.length = 0;
        /* initializing vector V */
        V.point[0].x = start.x;
        V.point[0].y = start.y;
        V.point[0].index = start.index;
        V.i = 0;
        V.j = 0;
        V.length = 0;
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
            V.point[1].x = points[i].x;
            V.point[1].y = points[i].y;
            V.point[1].index = points[i].index;
            V.i = V.point[1].x - V.point[0].x;
            V.j = V.point[1].y - V.point[0].y;
            V.length = length_v(V);
            /* initializing vector T2 */
            T2.point[1].x = V.point[1].x;
            T2.point[1].y = V.point[1].y;
            T2.point[1].index = INT_MAX;
            T2.i = (T2.point[1].x - T2.point[0].x) / V.length;
            T2.j = (T2.point[1].y - T2.point[0].y) / V.length;
            T2.point[1].x = V.point[0].x + T2.i;
            T2.point[1].y = V.point[0].y + T2.j;
            T2.length = length_v(T2);
            /* initializing tao, theta, and curvature */
            curr[i].tao = (dot_product(T1, T2)); //length of T1 and T2 is always 1
            if(curr[i].tao <= -1.0) {
                curr[i].tao = -1.0;
            }
            else if(curr[i].tao >= 1.0) {
                curr[i].tao = 1.0;
            }
            curr[i].x = V.point[1].x;
            curr[i].y = V.point[1].y;
            curr[i].index = V.point[1].index;
            curr[i].theta = angle_t(curr[i].tao);
            curr[i].curvature = calculate_curvature(T1, T2, curr[i].tao);
            curr[i].tao_distance = tao_distance(V, curr[i].curvature, curr[i].theta);
            V.point[1].tao_distance = curr[i].tao_distance;
            i++;
            count++;
        }
        /* sets the previous point as the previous best point */
        prev = best;
        /* find point with the lowest tao-distance */
        for(i = 0; i < size; i++) {
            if(best.tao_distance > curr[i].tao_distance) {
                best.x = curr[i].x;
                best.y = curr[i].y;
                best.index = curr[i].index;
                best.theta = curr[i].theta;
                best.curvature = curr[i].curvature;
                best.tao_distance = curr[i].tao_distance;
                k = i;
            }
        }
        /* record path */
        loop[m++] = k;
        /* if the best point has been visited before */
        if(visited[best.index] == 1) {
            m--;
            if(m != 0) {
                /* check if the loop has finished elsewhere
                if(loop[m] != loop[j]) {
                    for(j = m; j > 0; j--) {
                        if(loop[j] != loop[0]) {
                            m--;
                        }
                        else {
                            break;
                        }
                    }
                }*/
                for(j = 0; j < m; j++) {
                    printf("%d->", points[loop[j]].index);
                    mapped[loop[j]] = 1;
                }
                printf("%d\n", points[loop[m]].index);
                for(j = 0; j < m; j++) {
                    /* record segment */
                    tmp_segments[j][0] = loop[j];
                    tmp_segments[j][1] = loop[j + 1];
                }
                //tmp_segments[m][0] = loop[m];
                //tmp_segments[m][1] = loop[0];
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
                    printf("%d = (%d, %d): <%d,%d>\n", j, tmp_segments[j][0], tmp_segments[j][1], points[tmp_segments[j][0]].index, points[tmp_segments[j][1]].index);
                }
            }
            return;
        }
        /* reinitializing vector V */
        V.point[1].x = best.x;
        V.point[1].y = best.y;
        V.point[1].index = best.index;
        V.i = V.point[1].x - V.point[0].x;
        V.j = V.point[1].y - V.point[0].y;
        V.length = length_v(V);
        /* reinitializing vector T1 */
        T2.point[1].x = best.x;
        T2.point[1].y = best.y;
        T2.point[1].index = INT_MAX;
        T2.i = (T2.point[1].x - T2.point[0].x) / V.length;
        T2.j = (T2.point[1].y - T2.point[0].y) / V.length;
        T2.length = length_v(T2);
        T1.point[0].x = best.x;
        T1.point[0].y = best.y;
        T1.point[0].index = best.index;
        T1.point[1].x = best.x + T2.i;
        T1.point[1].y = best.y + T2.j;
        T1.point[1].index = INT_MAX;
        T1.i = (T1.point[1].x - T1.point[0].x);
        T1.j = (T1.point[1].y - T1.point[0].y);
        T1.length = length_v(T1);
        /* shifts starting point to best point */
        start.x = best.x;
        start.y = best.y;
        start.index = best.index;
        /* initializing vector T2 */
        T2.point[0].x = start.x;
        T2.point[0].y = start.y;
        T2.point[0].index = start.index;
        T2.i = 0;
        T2.j = 0;
        T2.length = 0;
        /* initializing vector V */
        V.point[0].x = start.x;
        V.point[0].y = start.y;
        V.point[0].index = start.index;
        V.i = 0;
        V.j = 0;
        V.length = 0;
        count = 0;
        permutations++;
        visited_count++;
    }
    return;
}

/* calculates the connections for un-joined vertices */
void join_vertex(vector<int *> *segments, struct point_t *points, struct point_t begin, int n, int size, FILE *gnu_files[NUM_FILES])
{
    struct vector_t V;
    struct vector_t T1;
    struct vector_t T2;
    struct point_t *curr = new struct point_t [size];
    struct point_t best;
    struct point_t start;
    struct point_t prev;
    struct point_t center;
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

    start.x = begin.x;
    start.y = begin.y;
    start.index = begin.index;
    prev.x = begin.x;
    prev.y = begin.y;
    prev.index = begin.index;
    best.x = begin.x;
    best.y = begin.y;
    best.index = begin.index;
    /* initializing vector T1 */
    T1.point[0].x = start.x;
    T1.point[0].y = start.y;
    T1.point[0].index = start.index;
    T1.i = (center.x - start.x) / distance_p(center, start);
    T1.j = (center.y - start.y) / distance_p(center, start);
    T1.point[1].x = start.x + T1.i;
    T1.point[1].y = start.y + T1.j;
    T1.point[1].index = INT_MAX;
    T1.length = length_v(T1);
    /* outer loop */
    while(m < 2) {
        /* store start index in visited-array */
        visited[start.index] = 1;
        i = 0;
        /* refreshing best index */
        best.tao_distance = DBL_MAX;
        best.index = start.index;
        /* initializing vector T2 */
        T2.point[0].x = start.x;
        T2.point[0].y = start.y;
        T2.point[0].index = start.index;
        T2.i = 0;
        T2.j = 0;
        T2.length = 0;
        /* initializing vector V */
        V.point[0].x = start.x;
        V.point[0].y = start.y;
        V.point[0].index = start.index;
        V.i = 0;
        V.j = 0;
        V.length = 0;
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
            V.point[1].x = points[i].x;
            V.point[1].y = points[i].y;
            V.point[1].index = points[i].index;
            V.i = V.point[1].x - V.point[0].x;
            V.j = V.point[1].y - V.point[0].y;
            V.length = length_v(V);
            /* initializing vector T2 */
            T2.point[1].x = V.point[1].x;
            T2.point[1].y = V.point[1].y;
            T2.point[1].index = INT_MAX;
            T2.i = (T2.point[1].x - T2.point[0].x) / V.length;
            T2.j = (T2.point[1].y - T2.point[0].y) / V.length;
            T2.point[1].x = V.point[0].x + T2.i;
            T2.point[1].y = V.point[0].y + T2.j;
            T2.length = length_v(T2);
            /* initializing tao, theta, and curvature */
            curr[i].tao = (dot_product(T1, T2)); //length of T1 and T2 is always 1
            if(curr[i].tao <= -1.0) {
                curr[i].tao = -1.0;
            }
            else if(curr[i].tao >= 1.0) {
                curr[i].tao = 1.0;
            }
            curr[i].x = V.point[1].x;
            curr[i].y = V.point[1].y;
            curr[i].index = V.point[1].index;
            curr[i].theta = angle_t(curr[i].tao);
            curr[i].curvature = calculate_curvature(T1, T2, curr[i].tao);
            curr[i].tao_distance = tao_distance(V, curr[i].curvature, curr[i].theta);
            V.point[1].tao_distance = curr[i].tao_distance;
            i++;
            count++;
        }
        /* sets the previous point as the previous best point */
        prev = best;
        /* find point with the lowest tao-distance */
        for(i = 0; i < size; i++) {
            if(best.tao_distance > curr[i].tao_distance) {
                best.x = curr[i].x;
                best.y = curr[i].y;
                best.index = curr[i].index;
                best.theta = curr[i].theta;
                best.curvature = curr[i].curvature;
                best.tao_distance = curr[i].tao_distance;
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
        /* reinitializing vector V */
        V.point[1].x = best.x;
        V.point[1].y = best.y;
        V.point[1].index = best.index;
        V.i = V.point[1].x - V.point[0].x;
        V.j = V.point[1].y - V.point[0].y;
        V.length = length_v(V);
        /* reinitializing vector T1 */
        T2.point[1].x = best.x;
        T2.point[1].y = best.y;
        T2.point[1].index = INT_MAX;
        T2.i = (T2.point[1].x - T2.point[0].x) / V.length;
        T2.j = (T2.point[1].y - T2.point[0].y) / V.length;
        T2.length = length_v(T2);
        T1.point[0].x = best.x;
        T1.point[0].y = best.y;
        T1.point[0].index = best.index;
        T1.point[1].x = best.x + T2.i;
        T1.point[1].y = best.y + T2.j;
        T1.point[1].index = INT_MAX;
        T1.i = (T1.point[1].x - T1.point[0].x);
        T1.j = (T1.point[1].y - T1.point[0].y);
        T1.length = length_v(T1);
        /* shifts starting point to best point */
        start.x = best.x;
        start.y = best.y;
        start.index = best.index;
        /* initializing vector T2 */
        T2.point[0].x = start.x;
        T2.point[0].y = start.y;
        T2.point[0].index = start.index;
        T2.i = 0;
        T2.j = 0;
        T2.length = 0;
        /* initializing vector V */
        V.point[0].x = start.x;
        V.point[0].y = start.y;
        V.point[0].index = start.index;
        V.i = 0;
        V.j = 0;
        V.length = 0;
        count = 0;
        permutations++;
    }
    return;
}

/* calculates the connections for un-joined segments */
void join_segment(vector<int *> *segments, struct point_t *points, struct point_t begin, struct point_t end, int n, int m, int size, FILE *gnu_files[NUM_FILES])
{
    struct vector_t V;
    struct vector_t T1;
    struct vector_t T2;
    struct point_t *curr = new struct point_t [size];
    struct point_t best;
    struct point_t start;
    struct point_t prev;
    struct point_t center;
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
    best.tao_distance = DBL_MAX;
    best.x = begin.x;
    best.y = begin.y;
    best.index = begin.index;
    /* initializing vector T1 */
    T1.point[0].x = begin.x;
    T1.point[0].y = begin.y;
    T1.point[0].index = begin.index;
    T1.i = (end.x - begin.x) / distance_p(end, begin);
    T1.j = (end.y - begin.y) / distance_p(end, begin);
    T1.point[1].x = begin.x + T1.i;
    T1.point[1].y = begin.y + T1.j;
    T1.point[1].index = INT_MAX;
    T1.length = length_v(T1);
    /* store start index in visited-array */
    visited[begin.index] = 1;
    i = 0;
    /* refreshing best index */
    best.tao_distance = DBL_MAX;
    best.index = begin.index;
    /* initializing vector T2 */
    T2.point[0].x = begin.x;
    T2.point[0].y = begin.y;
    T2.point[0].index = begin.index;
    T2.i = 0;
    T2.j = 0;
    T2.length = 0;
    /* initializing vector V */
    V.point[0].x = begin.x;
    V.point[0].y = begin.y;
    V.point[0].index = begin.index;
    V.i = 0;
    V.j = 0;
    V.length = 0;
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
        V.point[1].x = points[i].x;
        V.point[1].y = points[i].y;
        V.point[1].index = points[i].index;
        V.i = V.point[1].x - V.point[0].x;
        V.j = V.point[1].y - V.point[0].y;
        V.length = length_v(V);
        /* initializing vector T2 */
        T2.point[1].x = V.point[1].x;
        T2.point[1].y = V.point[1].y;
        T2.point[1].index = INT_MAX;
        T2.i = (T2.point[1].x - T2.point[0].x) / V.length;
        T2.j = (T2.point[1].y - T2.point[0].y) / V.length;
        T2.point[1].x = V.point[0].x + T2.i;
        T2.point[1].y = V.point[0].y + T2.j;
        T2.length = length_v(T2);
        /* initializing tao, theta, and curvature */
        curr[i].tao = (dot_product(T1, T2)); //length of T1 and T2 is always 1
        if(curr[i].tao <= -1.0) {
            curr[i].tao = -1.0;
        }
        else if(curr[i].tao >= 1.0) {
            curr[i].tao = 1.0;
        }
        curr[i].x = V.point[1].x;
        curr[i].y = V.point[1].y;
        curr[i].index = V.point[1].index;
        curr[i].theta = angle_t(curr[i].tao);
        curr[i].curvature = calculate_curvature(T1, T2, curr[i].tao);
        curr[i].tao_distance = tao_distance(V, curr[i].curvature, curr[i].theta) * curr[i].theta; //addition of multiply by theta
        V.point[1].tao_distance = curr[i].tao_distance;
        i++;
        count++;
    }
    /* find point with the lowest tao-distance */
    for(i = 0; i < size; i++) {
        if(best.tao_distance > curr[i].tao_distance) {
            best.x = curr[i].x;
            best.y = curr[i].y;
            best.index = curr[i].index;
            best.theta = curr[i].theta;
            best.curvature = curr[i].curvature;
            best.tao_distance = curr[i].tao_distance;
            k = i;
        }
    }
    /* record segment */
    tmp_segments[m][0] = m;
    tmp_segments[m][1] = k;
    printf("... joining: <%d,%d>\n", points[n].index, points[k].index);
    /* pushes new segments */
    pushed_segment = new int [2];
    pushed_segment[0] = tmp_segments[m][0];
    pushed_segment[1] = tmp_segments[m][1];
    segments->push_back(pushed_segment);
    return;
}

vector<struct polygon_t> construct_polygons(vector<int *> segments, struct point_t *points, int size, FILE *gnu_files[NUM_FILES])
{
    vector<vector<int> > tesselations;
    vector<struct polygon_t> polygons;
    struct polygon_t polygon;
    double sum_x = 0.0;
    double sum_y = 0.0;
    int i = 0;

    /* right-right addition */
    tesselations = tesselate(tesselations, segments, points, size, 'r', 'r', gnu_files);
    /* right-left addition */
    tesselations = tesselate(tesselations, segments, points, size, 'r', 'l', gnu_files);
    /* left-left addition */
    tesselations = tesselate(tesselations, segments, points, size, 'l', 'l', gnu_files);
    /* left-right addition */
    tesselations = tesselate(tesselations, segments, points, size, 'l', 'r', gnu_files);

    /* stores in polygon_t structure format */
    for(i = 0; i < tesselations.size(); i++) {
        polygon.shape = tesselations[i];
        polygon.perimeter = find_perimeter(tesselations[i], points);
        polygons.push_back(polygon);
    }

    return polygons;
}

/* calculate tesselations of polygons given all contoured segments */
vector<vector<int> > tesselate(vector<vector<int> > tesselations, vector<int *> segments, struct point_t *points, int size, char init, char add, FILE *gnu_files[NUM_FILES])
{
    struct vector_t X; //X-axis vector
    X.name = new char [2];
    X.name[0] = 'X';
    X.name[1] = '\0';
    struct vector_t Y; //Y-axis vector
    Y.name = new char [2];
    Y.name[0] = 'Y';
    Y.name[1] = '\0';
    struct point_t start;
    struct point_t prev;
    struct point_t root;
    vector<int *> edges; //pointer contains the index followed by position
    vector<int> path; //contains the path of nodes to add as a polygon
    vector<int> remove; //contains edges to remove
    int i = 0;
    int j = 0;
    int complete = 0;

    for(i = 0; i < size; i++) {
        point_t start = points[i];
        point_t root = points[i];
        /* find the initial cluster of edges */
        edges = edge_search(segments, start.index, points, size);
        /* initialize axis vectors */
        Y.point[0].x = start.x;
        Y.point[0].y = start.y;
        Y.point[0].index = start.index;
        Y.i = 0;
        Y.j = 1;
        Y.point[1].x = start.x + Y.i;
        Y.point[1].y = start.y + Y.j;
        Y.point[1].index = -1;
        Y.length = 1;
        X.point[0].x = start.x;
        X.point[0].y = start.y;
        X.point[0].index = start.index;
        X.i = 1;
        X.j = 0;
        X.point[1].x = start.x + X.i;
        X.point[1].y = start.y + X.j;
        X.point[1].index = -1;
        X.length = 1;
        /* prints the initial cluster of edges */
        printf("\n");
        printf("TYPE: %c %c | ", init, add);
        printf("INITIAL: ");
        for(j = 0; j < edges.size(); j++) {
            printf("%d ", points[edges[j][1]].index);
        }
        printf("| ");
        /* find the initial direction */
        path.push_back(edges[0][0]);
        path = init_path(path, edges, points, X, Y, init);
        if(path.size() == 0) {
            continue;
        }
        printf("PATH: %d %d ", points[path[0]].index, points[path[1]].index);
        j = 1;
        complete = 0;
        while(!complete) {
            prev = start;
            start = points[path[j]];
            /* initialize axis vectors */
            Y.point[0].x = start.x;
            Y.point[0].y = start.y;
            Y.point[0].index = start.index;
            Y.i = (start.x - prev.x) / distance_p(prev, start);
            Y.j = (start.y - prev.y) / distance_p(prev, start);
            Y.point[1].x = start.x + Y.i;
            Y.point[1].y = start.y + Y.j;
            Y.point[1].index = -1;
            Y.length = length_v(Y);
            X.point[0].x = start.x;
            X.point[0].y = start.y;
            X.point[0].index = start.index;
            X.i = Y.j;
            X.j = -Y.i;
            X.point[1].x = start.x + X.i;
            X.point[1].y = start.y + X.j;
            X.point[1].index = -1;
            X.length = length_v(X);
            /* find the initial tesselations of edges */
            edges = edge_search(segments, start.index, points, size);
            path = add_path(path, edges, points, X, Y, add);
            if(path.size() == 0) {
                break;
            }
            printf("%d ", points[path[j + 1]].index);
            if(path[j + 1] == path[0]) {
                complete = 1;
            }
            j++;
        }
        if(!complete) {
            path.clear();
            remove.clear();
            continue;
        }
        printf("\n\n");
        tesselations.push_back(path);
        /* find the initial tesselations of edges */
        edges = edge_search(segments, root.index, points, size);
        /* find the index in edges that matches path[1] */
        if(index_match(edges, path[1]) > -1) {
            remove.push_back(index_match(edges, path[1]));
        }
        sort(remove.begin(), remove.end());
        /* erasing the edges already traversed */
        for(j = 0; j < remove.size(); j++) {
            edges.erase(edges.begin() + remove[j] - j);
        }
        /* erasing the ending edge, no backtracking will occur */
        if(index_match(edges, path[path.size() - 2]) > -1) {
            edges.erase(edges.begin() + index_match(edges, path[path.size() - 2]));
        }
        path.clear();
        remove.clear();
    }
    return tesselations;
}

/* initializes the path from a vertex */
vector<int> init_path(vector<int> path, vector<int *> edges, struct point_t *points, struct vector_t X, struct vector_t Y, char type)
{
    struct vector_t E; //edge vector
    E.name = new char [2];
    E.name[0] = 'E';
    E.name[1] = '\0';
    double *X_flags = new double [edges.size()];
    double *Y_flags = new double [edges.size()];
    double dot = 0.0;
    int *quads = new int [edges.size()];
    int i = 0;
    int curr = -1;
    int init = 0;
    int found = 0;
    /* calculate and set flags for each edge */
    for(i = 0; i < edges.size(); i++) {
        /* set edge vector */
        E.point[0].x = points[edges[i][0]].x;
        E.point[0].y = points[edges[i][0]].y;
        E.point[0].index = points[edges[i][0]].index;
        E.point[1].x = points[edges[i][1]].x;
        E.point[1].y = points[edges[i][1]].y;
        E.point[1].index = points[edges[i][1]].index;
        //E.i = (points[edges[i][1]].x - points[edges[i][0]].x) / distance_p(points[edges[i][1]], points[edges[i][0]]);
        //E.j = (points[edges[i][1]].y - points[edges[i][0]].y) / distance_p(points[edges[i][1]], points[edges[i][0]]);
        E.i = (points[edges[i][1]].x - points[edges[i][0]].x);
        E.j = (points[edges[i][1]].y - points[edges[i][0]].y);
        E.length = length_v(E);
        X_flags[i] = dot_product(E, X);
        Y_flags[i] = dot_product(E, Y);
        printf("\n");
        print_v(E);
        print_v(X);
        print_v(Y);
        printf("DOTX: %0.2lf, ", X_flags[i]);
        printf("DOTY: %0.2lf, ", Y_flags[i]);
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
            else {
                quads[i] = 6;
            }
        }
        printf("QUAD: %d\n", quads[i]);
    }
    /* check for which direction to initialize to */
    switch(type) {
    case 'r':
        /* first check for an edge in both quadrants 1 and 2 */
        for(i = 0; i < edges.size(); i++) {
            if(quads[i] == 5) {
                curr = i;
                found = 1;
                break;
            }
        }
        /* check for an edge in quadrant 1 */
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
        /* check for an edge in quadrant 4 */
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
        /* check for an edge in both quadrants 3 and 4 */
        if(!found) {
            for(i = 0; i < edges.size(); i++) {
                if(quads[i] == 6) {
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
                curr = i;
                found = 1;
                break;
            }
        }
        /* first check for an edge in quadrant 2 */
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
        /* now check for an edge in quadrant 3 */
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
        /* now check for an edge in both quadrants 3 and 4 */
        if(!found) {
            for(i = 0; i < edges.size(); i++) {
                if(quads[i] == 6) {
                    curr = i;
                    break;
                }
            }
        }
        break;
    default:
        printf("\nDIRECTION NOT SET!\n\n");
        exit(EXIT_FAILURE);
    }
    /* use curr to create the link in the path */
    if(curr > -1) {
        path.push_back(edges[curr][1]);
        if(duplicate_search(path)) {
            printf("duplicate found: %d\n", edges[curr][1]);
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
vector<int> add_path(vector<int> path, vector<int *> edges, struct point_t *points, struct vector_t X, struct vector_t Y, char type)
{
    struct vector_t E; //edge vector
    E.name = new char [2];
    E.name[0] = 'E';
    E.name[1] = '\0';
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
        /* set edge vector */
        E.point[0].x = points[edges[i][0]].x;
        E.point[0].y = points[edges[i][0]].y;
        E.point[0].index = points[edges[i][0]].index;
        E.point[1].x = points[edges[i][1]].x;
        E.point[1].y = points[edges[i][1]].y;
        E.point[1].index = points[edges[i][1]].index;
        E.i = (points[edges[i][1]].x - points[edges[i][0]].x) / distance_p(points[edges[i][1]], points[edges[i][0]]);
        E.j = (points[edges[i][1]].y - points[edges[i][0]].y) / distance_p(points[edges[i][1]], points[edges[i][0]]);
        E.length = length_v(E);
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
            /* if the segement is directly behind then stop */
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
        printf("\nDIRECTION NOT SET!\n\n");
        exit(EXIT_FAILURE);
    }
    /* use curr to create the link in the path */
    if(curr > -1) {
        path.push_back(edges[curr][1]);
        /* if path is not a viable loop */
        if(duplicate_search(path)) {
            printf("duplicate found: %d\n", points[edges[curr][1]].index);
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
vector<int *> edge_search(vector<int *> segments, int vertex, struct point_t *points, int size)
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

/* return the index of the matching element */
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

/* searches through a shape for a matching vertex */
int shape_search(vector<int> shape, int vertex)
{
    /* check if the vertex is found */
    if(find(shape.begin(), shape.end(), vertex) != shape.end()) {
        return 1;
    }
    else {
        return 0;
    }
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
double find_perimeter(vector<int> shape, struct point_t *points)
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
        printf("new: ");
        for(j = 0; j < (tmp[i]).shape.size(); j++) {
            printf("%d ", tmp[i].shape[j]);
        }
        printf("\n");
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
                    printf("%d. deleting: ", j);
                    for(k = 0; k < (tmp[j]).shape.size(); k++) {
                        printf("%d ", tmp[j].shape[k]);
                    }
                    printf("\n");
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
vector<struct polygon_t> optimize_polygons(vector<struct polygon_t> polygons, vector<int *> *segments, struct point_t *points, int size, FILE *gnu_files[NUM_FILES])
{
    struct vector_t V;
    struct vector_t T1;
    struct vector_t T2;
    struct point_t *curr = new struct point_t [size];
    struct point_t best;
    struct point_t start;
    struct point_t end;
    struct point_t prev;
    struct point_t center;
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
        start.x = points[(polygons[i]).shape[0]].x;
        start.y = points[(polygons[i]).shape[0]].y;
        start.index = points[(polygons[i]).shape[0]].index;
        prev.x = start.x;
        prev.y = start.y;
        prev.index = start.index;
        best.x = start.x;
        best.y = start.y;
        best.index = start.index;
        /* initializing vector T1 */
        T1.point[0].x = start.x;
        T1.point[0].y = start.y;
        T1.point[0].index = start.index;
        T1.i = (center.x - start.x) / distance_p(center, start);
        T1.j = (center.y - start.y) / distance_p(center, start);
        T1.point[1].x = start.x + T1.i;
        T1.point[1].y = start.y + T1.j;
        T1.point[1].index = INT_MAX;
        T1.length = length_v(T1);
        /* loops through all segments of the polygon */
        for(j = 0; j < (polygons[i]).shape.size() - 2; j++) {
            k = (polygons[i]).shape[j];
            l = (polygons[i]).shape[j + 1];
            start.x = points[k].x;
            start.y = points[k].y;
            start.index = points[k].index;
            end.x = points[l].x;
            end.y = points[l].y;
            end.index = points[l].index;
            /* initializing vector T1 */
            T1.point[0].x = start.x;
            T1.point[0].y = start.y;
            T1.point[0].index = start.index;
            T1.i = (end.x - start.x) / distance_p(end, start);
            T1.j = (end.y - start.y) / distance_p(end, start);
            T1.point[1].x = start.x + T1.i;
            T1.point[1].y = start.y + T1.j;
            T1.point[1].index = INT_MAX;
            T1.length = length_v(T1);
            /* refreshing best index */
            best.tao_distance = DBL_MAX;
            /* initializing vector T2 */
            T2.point[0].x = start.x;
            T2.point[0].y = start.y;
            T2.point[0].index = start.index;
            T2.i = 0;
            T2.j = 0;
            T2.length = 0;
            /* initializing vector V */
            V.point[0].x = start.x;
            V.point[0].y = start.y;
            V.point[0].index = start.index;
            V.i = 0;
            V.j = 0;
            V.length = 0;
            /* loops through all possible indices from start */
            for(k = 0; k < size; k++) {
                /* initializing vector V */
                V.point[1].x = points[k].x;
                V.point[1].y = points[k].y;
                V.point[1].index = points[k].index;
                V.i = V.point[1].x - V.point[0].x;
                V.j = V.point[1].y - V.point[0].y;
                V.length = length_v(V);
                /* initializing vector T2 */
                T2.point[1].x = V.point[1].x;
                T2.point[1].y = V.point[1].y;
                T2.point[1].index = INT_MAX;
                T2.i = (T2.point[1].x - T2.point[0].x) / V.length;
                T2.j = (T2.point[1].y - T2.point[0].y) / V.length;
                T2.point[1].x = V.point[0].x + T2.i;
                T2.point[1].y = V.point[0].y + T2.j;
                T2.length = length_v(T2);
                /* initializing tao, theta, and curvature */
                curr[k].tao = (dot_product(T1, T2)); //length of T1 and T2 is always 1
                if(curr[k].tao <= -1.0) {
                    curr[k].tao = -1.0;
                }
                else if(curr[k].tao >= 1.0) {
                    curr[k].tao = 1.0;
                }
                curr[k].x = V.point[1].x;
                curr[k].y = V.point[1].y;
                curr[k].index = V.point[1].index;
                curr[k].theta = angle_t(curr[k].tao);
                curr[k].curvature = calculate_curvature(T1, T2, curr[k].tao);
                curr[k].tao_distance = tao_distance(V, curr[k].curvature, curr[k].theta);
                V.point[1].tao_distance = curr[k].tao_distance;
            }
            /* find point with the lowest tao-distance */
            for(k = 0; k < size; k++) {
                if(best.tao_distance > curr[k].tao_distance) {
                    best.x = curr[k].x;
                    best.y = curr[k].y;
                    best.index = curr[k].index;
                    best.theta = curr[k].theta;
                    best.curvature = curr[k].curvature;
                    best.tao_distance = curr[k].tao_distance;
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
                    printf("added: <%d,%d>\n", points[(polygons[i]).shape[j]].index, points[n].index);
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
        start.x = points[(polygons[i]).shape[0]].x;
        start.y = points[(polygons[i]).shape[0]].y;
        start.index = points[(polygons[i]).shape[0]].index;
        prev.x = start.x;
        prev.y = start.y;
        prev.index = start.index;
        best.x = start.x;
        best.y = start.y;
        best.index = start.index;
        /* initializing vector T1 */
        T1.point[0].x = start.x;
        T1.point[0].y = start.y;
        T1.point[0].index = start.index;
        T1.i = (center.x - start.x) / distance_p(center, start);
        T1.j = (center.y - start.y) / distance_p(center, start);
        T1.point[1].x = start.x + T1.i;
        T1.point[1].y = start.y + T1.j;
        T1.point[1].index = INT_MAX;
        T1.length = length_v(T1);
        /* loops through all indices of a single polygon */
        for(j = 0; j < (polygons[i]).shape.size(); j++) {
            k = (polygons[i]).shape[j];
            start.x = points[k].x;
            start.y = points[k].y;
            start.index = points[k].index;
            /* initializing vector T1 */
            T1.point[0].x = start.x;
            T1.point[0].y = start.y;
            T1.point[0].index = start.index;
            T1.i = (center.x - start.x) / distance_p(center, start);
            T1.j = (center.y - start.y) / distance_p(center, start);
            T1.point[1].x = start.x + T1.i;
            T1.point[1].y = start.y + T1.j;
            T1.point[1].index = INT_MAX;
            T1.length = length_v(T1);
            /* refreshing best index */
            best.tao_distance = DBL_MAX;
            /* initializing vector T2 */
            T2.point[0].x = start.x;
            T2.point[0].y = start.y;
            T2.point[0].index = start.index;
            T2.i = 0;
            T2.j = 0;
            T2.length = 0;
            /* initializing vector V */
            V.point[0].x = start.x;
            V.point[0].y = start.y;
            V.point[0].index = start.index;
            V.i = 0;
            V.j = 0;
            V.length = 0;
            /* loops through all possible indices from start */
            for(k = 0; k < size; k++) {
                /* initializing vector V */
                V.point[1].x = points[k].x;
                V.point[1].y = points[k].y;
                V.point[1].index = points[k].index;
                V.i = V.point[1].x - V.point[0].x;
                V.j = V.point[1].y - V.point[0].y;
                V.length = length_v(V);
                /* initializing vector T2 */
                T2.point[1].x = V.point[1].x;
                T2.point[1].y = V.point[1].y;
                T2.point[1].index = INT_MAX;
                T2.i = (T2.point[1].x - T2.point[0].x) / V.length;
                T2.j = (T2.point[1].y - T2.point[0].y) / V.length;
                T2.point[1].x = V.point[0].x + T2.i;
                T2.point[1].y = V.point[0].y + T2.j;
                T2.length = length_v(T2);
                /* initializing tao, theta, and curvature */
                curr[k].tao = (dot_product(T1, T2)); //length of T1 and T2 is always 1
                if(curr[k].tao <= -1.0) {
                    curr[k].tao = -1.0;
                }
                else if(curr[k].tao >= 1.0) {
                    curr[k].tao = 1.0;
                }
                curr[k].x = V.point[1].x;
                curr[k].y = V.point[1].y;
                curr[k].index = V.point[1].index;
                curr[k].theta = angle_t(curr[k].tao);
                curr[k].curvature = calculate_curvature(T1, T2, curr[k].tao);
                curr[k].tao_distance = tao_distance(V, curr[k].curvature, curr[k].theta);
                V.point[1].tao_distance = curr[k].tao_distance;
            }
            /* find point with the lowest tao-distance */
            for(k = 0; k < size; k++) {
                if(best.tao_distance > curr[k].tao_distance) {
                    best.x = curr[k].x;
                    best.y = curr[k].y;
                    best.index = curr[k].index;
                    best.theta = curr[k].theta;
                    best.curvature = curr[k].curvature;
                    best.tao_distance = curr[k].tao_distance;
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
                    printf("added: <%d,%d>\n", points[(polygons[i]).shape[j]].index, points[n].index);
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
        start.x = points[(polygons[i]).shape[0]].x;
        start.y = points[(polygons[i]).shape[0]].y;
        start.index = points[(polygons[i]).shape[0]].index;
        prev.x = start.x;
        prev.y = start.y;
        prev.index = start.index;
        best.x = start.x;
        best.y = start.y;
        best.index = start.index;
        /* initializing vector T1 */
        T1.point[0].x = start.x;
        T1.point[0].y = start.y;
        T1.point[0].index = start.index;
        T1.i = -(center.x - start.x) / distance_p(center, start);
        T1.j = -(center.y - start.y) / distance_p(center, start);
        T1.point[1].x = start.x + T1.i;
        T1.point[1].y = start.y + T1.j;
        T1.point[1].index = INT_MAX;
        T1.length = length_v(T1);
        /* loops through all indices of a single polygon */
        for(j = 0; j < (polygons[i]).shape.size(); j++) {
            k = (polygons[i]).shape[j];
            start.x = points[k].x;
            start.y = points[k].y;
            start.index = points[k].index;
            /* initializing vector T1 */
            T1.point[0].x = start.x;
            T1.point[0].y = start.y;
            T1.point[0].index = start.index;
            T1.i = -(center.x - start.x) / distance_p(center, start);
            T1.j = -(center.y - start.y) / distance_p(center, start);
            T1.point[1].x = start.x + T1.i;
            T1.point[1].y = start.y + T1.j;
            T1.point[1].index = INT_MAX;
            T1.length = length_v(T1);
            /* refreshing best index */
            best.tao_distance = DBL_MAX;
            /* initializing vector T2 */
            T2.point[0].x = start.x;
            T2.point[0].y = start.y;
            T2.point[0].index = start.index;
            T2.i = 0;
            T2.j = 0;
            T2.length = 0;
            /* initializing vector V */
            V.point[0].x = start.x;
            V.point[0].y = start.y;
            V.point[0].index = start.index;
            V.i = 0;
            V.j = 0;
            V.length = 0;
            /* loops through all possible indices from start */
            for(k = 0; k < size; k++) {
                /* initializing vector V */
                V.point[1].x = points[k].x;
                V.point[1].y = points[k].y;
                V.point[1].index = points[k].index;
                V.i = V.point[1].x - V.point[0].x;
                V.j = V.point[1].y - V.point[0].y;
                V.length = length_v(V);
                /* initializing vector T2 */
                T2.point[1].x = V.point[1].x;
                T2.point[1].y = V.point[1].y;
                T2.point[1].index = INT_MAX;
                T2.i = (T2.point[1].x - T2.point[0].x) / V.length;
                T2.j = (T2.point[1].y - T2.point[0].y) / V.length;
                T2.point[1].x = V.point[0].x + T2.i;
                T2.point[1].y = V.point[0].y + T2.j;
                T2.length = length_v(T2);
                /* initializing tao, theta, and curvature */
                curr[k].tao = (dot_product(T1, T2)); //length of T1 and T2 is always 1
                if(curr[k].tao <= -1.0) {
                    curr[k].tao = -1.0;
                }
                else if(curr[k].tao >= 1.0) {
                    curr[k].tao = 1.0;
                }
                curr[k].x = V.point[1].x;
                curr[k].y = V.point[1].y;
                curr[k].index = V.point[1].index;
                curr[k].theta = angle_t(curr[k].tao);
                curr[k].curvature = calculate_curvature(T1, T2, curr[k].tao);
                curr[k].tao_distance = tao_distance(V, curr[k].curvature, curr[k].theta);
                V.point[1].tao_distance = curr[k].tao_distance;
            }
            /* find point with the lowest tao-distance */
            for(k = 0; k < size; k++) {
                if(best.tao_distance > curr[k].tao_distance) {
                    best.x = curr[k].x;
                    best.y = curr[k].y;
                    best.index = curr[k].index;
                    best.theta = curr[k].theta;
                    best.curvature = curr[k].curvature;
                    best.tao_distance = curr[k].tao_distance;
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
                    printf("added: <%d,%d>\n", points[(polygons[i]).shape[j]].index, points[n].index);
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
    /* make sure each 2-edge point is optimal
    for(i = 0; i < size; i++) {
        edges = edge_search(*segments, points[i].index, points, size);
        printf("EDGES: %d\n", edges.size());
        if(edges.size() == 2) {
            /* check that the point is optimal
            best.tao_distance = DBL_MAX;
            for(j = 0; j < size; j++) {
                if(i == j) {
                    continue;
                }
                curr->tao_distance = distance_p(points[i], points[j]);
                if(curr->tao_distance < best.tao_distance) {
                    best.tao_distance = curr->tao_distance;
                    n = j;
                }
            }
            /* record path
            for(k = 0; k < segments->size(); k++) {
                if(segment_match(*segments, points[i].index, points[n].index) == -1) {
                    tmp = new int [2];
                    tmp[0] = points[i].index;
                    tmp[1] = points[n].index;
                    segments->push_back(tmp);
                    printf("added: <%d,%d>\n", points[i].index, points[n].index);
                    break;
                }
            }
        }
    }*/
    return polygons;
}

/* returns the polygons with all crosses removed */
void remove_crosses(vector<int *> *segments, struct point_t *points, int size, FILE *gnu_files[NUM_FILES])
{
    struct point_t p1;
    struct point_t p2;
    struct point_t p3;
    struct point_t p4;
    struct point_t tmp;
    int i = 0;
    int j = 0;
    double y = 0.0;
    double x = 0.0;
    double m1 = 0.0;
    double m2 = 0.0;
    double b1 = 0.0;
    double b2 = 0.0;
    double r = DBL_MAX;
    double l = DBL_MIN;

    /* loops through all segments */
    for(i = 0; i < segments->size(); i++) {
        /* loops through all segments */
        for(j = 0; j < segments->size(); j++) {
            if(i == j) {
                continue;
            }
            printf("checking (%d,%d): <%d,%d> ", i, j, points[(*segments)[i][0]].index, points[(*segments)[i][1]].index);
            printf("and <%d,%d>\n", points[(*segments)[j][0]].index, points[(*segments)[j][1]].index);
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
            if(x < r && x > l) {
                printf("\nINTERSECTION FOUND!!!\n");
                /* remove the intersection */
                printf("intersection: (%0.2lf, %0.2lf), ", x, y);
                printf("<%d,%d> ", points[(*segments)[i][0]].index, points[(*segments)[i][1]].index);
                printf("<%d,%d>\n", points[(*segments)[j][0]].index, points[(*segments)[j][1]].index);
                if(i < j) {
                    segments->erase(segments->begin() + j);
                    segments->erase(segments->begin() + i);
                }
                else {
                    segments->erase(segments->begin() + i);
                    segments->erase(segments->begin() + j);
                }
                i = 0;
                j = 0;
            }
        }
    }
}

/* returns the index of the requested vertex */
int point_match(struct point_t *points, int size, int vertex)
{
    for(int i = 0; i < size; i++) {
       if(points[i].index == vertex) {
            return i;
        }
    }
    return -1;
}

/* calculates curvature given structure k */
double calculate_curvature(struct vector_t T1, struct vector_t T2, double tao)
{
    return (distance_v(T1, T2) / angle_t(tao));
}

/* calculates angle given tao */
double angle_t(double tao)
{
    return (acos(tao) + (M_PI / 180));
}

/* calculates angle between two vectors */
double angle_v(struct vector_t V1, struct vector_t V2)
{
    return (acos(dot_product(V1, V2) / (V1.length * V2.length)));
}

/* calculates distance given index and structure */
double tao_distance(struct vector_t V, double curvature, double theta)
{
    return (V.length + curvature + theta);
}

/* calculates distance given two points */
double distance_p(struct point_t start, struct point_t end)
{
    return sqrt(pow(end.x - start.x, 2) + pow(end.y - start.y, 2));
}

/* calculates distance given two vectors */
double distance_v(struct vector_t V1, struct vector_t V2)
{
    return sqrt(pow(V2.i - V1.i, 2) + pow(V2.j - V1.j, 2));
}

/* calculates length of a single vectors */
double length_v(struct vector_t V)
{
    return sqrt(pow(V.i, 2) + pow(V.j, 2));
}

/* calculates dot product of two vectors */
double dot_product(struct vector_t V1, struct vector_t V2)
{
    return ((V2.i * V1.i) + (V2.j * V1.j));
}

/* prints vector structure for debugging */
void print_v(struct vector_t V)
{
    printf("%s: datapoints = %d (%0.3lf,%0.3lf), %d (%0.3lf,%0.3lf)\n   components = <%0.3lf,%0.3lf>\n   length     = %0.3lf\n\n", V.name, V.point[0].index, V.point[0].x, V.point[0].y, V.point[1].index, V.point[1].x, V.point[1].y, V.i, V.j, V.length);
}

/* prints curvature structure for debugging */
void print(struct vector_t V, struct vector_t T1, struct vector_t T2, double curvature, double theta, double tao, double tao_distance)
{
    printf("V: %d[0](%lf, %lf), %d[1](%lf, %lf), <%lf, %lf>, |V| = %lf\n", V.point[0].index, V.point[0].x, V.point[0].y, V.point[1].index, V.point[1].x, V.point[1].y, V.i, V.j, V.length);
    printf("T1: point[0](%lf, %lf), point[1](%lf, %lf), <%lf, %lf>, |T1| = %lf\n", T1.point[0].x, T1.point[0].y, T1.point[1].x, T1.point[1].y, T1.i, T1.j, T1.length);
    printf("T2: point[0](%lf, %lf), point[1](%lf, %lf), <%lf, %lf>, |T2| = %lf\n", T2.point[0].x, T2.point[0].y, T2.point[1].x, T2.point[1].y, T2.i, T2.j, T2.length);
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
