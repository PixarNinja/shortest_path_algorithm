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

#define NUM_FILES 8

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

void construct_segments(vector<int *> *segments, struct point_t *points, struct point_t begin, int n, int size, FILE *gnu_files[NUM_FILES], int *mapped, int **recorded);
void join_vertex(vector<int *> *segments, struct point_t *points, struct point_t begin, int n, int size);
void join_segment(vector<int *> *segments, struct point_t *points, struct point_t begin, struct point_t end, int n, int m, int size);
vector<struct polygon_t> construct_polygons(vector<int *> segments, struct point_t *points, int size);
vector<vector<int> > tessellate(vector<vector<int> > tessellations, vector<int *> segments, struct point_t *points, int size, char init, char add, int branch);
vector<struct polygon_t> tessellate_cross(vector<int *> segments, int i, int j, struct point_t *points, int size);
vector<int> *find_shape(vector<int *> segments, struct point_t *points, struct point_t start, int size, char init, char add, struct vector_t X, struct vector_t Y);
vector<int> init_path(vector<int> path, vector<int *> edges, struct point_t *points, struct vector_t X, struct vector_t Y, char type, int branch);
vector<int> add_path(vector<int> path, vector<int *> edges, struct point_t *points, struct vector_t X, struct vector_t Y, char type);
vector<int *> edge_search(vector<int *> segments, int vertex, struct point_t *points, int size);
int index_match(vector<int *> segments, int vertex);
int shape_search(vector<int> shape, int vertex);
int edge_match(struct polygon_t polygon, int *edge);
int polygons_search(vector<vector<int> > polygons, int vertex);
double find_perimeter(vector<int> shape, struct point_t *points);
int segment_match(vector<int *> segments, int beginning, int end);
int duplicate_search(vector<int> shape);
vector<struct polygon_t> delete_duplicates(vector<struct polygon_t> polygons);
vector<struct polygon_t> optimize_polygons(vector<struct polygon_t> polygons, vector<int *> *segments, struct point_t *points, int size);
void remove_crosses(vector<int *> *segments, struct point_t *points, int size);
void finalize_segments(vector<int *> *segments, struct point_t *points, int size);
int intersection(struct point_t p1, struct point_t p2, struct point_t p3, struct point_t p4);
struct polygon_t find_shortest_path(vector<struct polygon_t> polygons, struct point_t *points, int size);
int accept_polygon(struct polygon_t polygon, vector<int *> segments, struct point_t *points);
int smallest_neighbour(vector<struct polygon_t> polygons, struct polygon_t source, int n);
vector<int *> disjoint_edges(struct polygon_t A, struct polygon_t B);
vector<int *> shared_edges(struct polygon_t A, struct polygon_t B);
vector<int> shared_points(struct polygon_t A, struct polygon_t B);
void visit_polygon(int *visited, struct polygon_t polygon, struct point_t *points);
struct polygon_t add_polygons(struct polygon_t A, struct polygon_t B, struct point_t *points);
struct polygon_t sub_polygons(struct polygon_t A, struct polygon_t B, struct point_t *points);
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
    if(argc != 3) {
        printf("\nUsage: ./tessellate [datapoint path] [output path]\n\n");
        exit(EXIT_FAILURE);
    }
    FILE *data;
    FILE *output;
    FILE *gnu_files[NUM_FILES];
    vector<struct polygon_t> polygons;
    struct polygon_t tmp_polygon;
    struct polygon_t shortest_path;
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
    /* grabs name for output files */
    char *name = new char [1024];
    char *tmp = new char [strlen(argv[1]) + 1];
    for(i = 0; i < strlen(argv[1]) + 1; i++)
        tmp[i] = argv[1][i];
    tmp = strtok(tmp, "/");
    while (tmp != NULL) {
        for(i = 0; i < strlen(tmp) + 1; i++)
            name[i] = tmp[i];
        tmp = strtok (NULL, "/");
    }
    /* setup pathnames for plot output */
    name = strtok(name, ".");
    const char *gnu_path = "./gnu_files/";
    char *commands = new char [strlen(gnu_path) + strlen(name) + strlen("_commands.gpf") + 1];
    char *datapoints = new char [strlen(gnu_path) + strlen(name) + strlen("_datapoints.gpf") + 1];
    char *lines = new char [strlen(gnu_path) + strlen(name) + strlen("_lines.gpf") + 1];
    char *extrapoints = new char [strlen(gnu_path) + strlen(name) + strlen("_lines.gpf") + 1];
    char *centerpoint = new char [strlen(gnu_path) + strlen(name) + strlen("_centerpoint.gpf") + 1];
    char *path = new char [strlen(gnu_path) + strlen(name) + strlen("_calculated_path.gpf") + 1];
    char *shortest = new char [strlen(gnu_path) + strlen(name) + strlen("_shortest_path.gpf") + 1];
    char *gnu_tmp = new char [strlen(gnu_path) + strlen("tmp.gpf") + 1];
    sprintf(commands, "%s%s_commands.gpf", gnu_path, name);
    sprintf(datapoints, "%s%s_datapoints.gpf", gnu_path, name);
    sprintf(lines, "%s%s_lines.gpf", gnu_path, name);
    sprintf(extrapoints, "%s%s_extrapoints.gpf", gnu_path, name);
    sprintf(centerpoint, "%s%s_centerpoint.gpf", gnu_path, name);
    sprintf(path, "%s%s_calculated_path.gpf", gnu_path, name);
    sprintf(shortest, "%s%s_shortest_path.gpf", gnu_path, name);
    sprintf(gnu_tmp, "%stmp.gpf", gnu_path);
    gnu_files[0] = fopen(commands, "w+");
    gnu_files[1] = fopen(datapoints, "w+");
    gnu_files[2] = fopen(lines, "w+");
    gnu_files[3] = fopen(extrapoints, "w+");
    gnu_files[4] = fopen(centerpoint, "w+");
    gnu_files[5] = fopen(path, "w+");
    gnu_files[6] = fopen(shortest, "r");
    gnu_files[7] = fopen(gnu_tmp, "w+");
    data = fopen(argv[1], "r");
    output = fopen(argv[2], "w+");
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
    /* plot datapoints */
    for(i = 0; i < size; i++) {
        fprintf(gnu_files[1], "%lf %lf %d\n", points[i].x, points[i].y, points[i].index);
    }
    /* plot setup */
    fprintf(gnu_files[0], "set xrange [%lf:%lf]\n", -(range + 1), range + 1);
    fprintf(gnu_files[0], "set yrange [%lf:%lf]\n", -(range + 1), range + 1);
    fprintf(gnu_files[0], "set size ratio 1\n");
    fprintf(gnu_files[0], "set grid\n");
    fprintf(gnu_files[0], "set title \"%s\"\n", argv[1]);
    fprintf(gnu_files[0], "set style line 1 lc rgb \"black\" lw 1\n");
    fprintf(gnu_files[0], "set style line 2 lc rgb \"red\" lw 3\n");
    fprintf(gnu_files[0], "set style line 3 lc rgb \"blue\" lw 2\n");
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
    /* get rid of crossing lines */
    remove_crosses(segments, points, size);
    /* optimize tessellations */
    do {
        /* find polygons */
        polygons = construct_polygons(*segments, points, size);
        polygons = delete_duplicates(polygons);
        /* ensure each polygon is optimal */
        count = segments->size();
        polygons = optimize_polygons(polygons, segments, points, size);
        /* get rid of crossing lines */
        remove_crosses(segments, points, size);
    } while(count < segments->size());
    /* finalize segments */
    finalize_segments(segments, points, size);
    /* find polygons */
    polygons = construct_polygons(*segments, points, size);
    polygons = delete_duplicates(polygons);
    /* bubble sort polygons by perimeter */
    for(i = 0; i < polygons.size(); i++) {
        for(j = polygons.size() - 1; j > i; j--) {
            if(polygons[j].perimeter < polygons[j - 1].perimeter) {
                tmp_polygon = polygons[j];
                polygons[j] = polygons[j - 1];
                polygons[j - 1] = tmp_polygon;
            }
        }
    }
    /* pop the incorrect polygon */
    polygons.pop_back();
    /* find shortest path
    shortest_path = find_shortest_path(polygons, points, size);
    /* plot shortest path
    printf("\nCALCULATED PATH: ");
    for(i = 0; i < shortest_path.shape.size(); i++) {
        printf("%d->", points[shortest_path.shape[i]].index);
        fprintf(gnu_files[5], "%lf %lf\n", points[shortest_path.shape[i]].x, points[shortest_path.shape[i]].y);
    }
    printf("%d\n\n", points[shortest_path.shape[0]].index);
    /* plot segment information */
    for(i = 0; i < segments->size(); i++) {
        fprintf(gnu_files[2], "%lf %lf %d\n", points[(*segments)[i][0]].x, points[(*segments)[i][0]].y, points[(*segments)[i][0]].index);
        fprintf(gnu_files[2], "%lf %lf %d\n", points[(*segments)[i][1]].x, points[(*segments)[i][1]].y, points[(*segments)[i][1]].index);
        fprintf(gnu_files[2], "\n");
        /* write segments to output file */
        fprintf(output, "%d %d\n", points[(*segments)[i][0]].index, points[(*segments)[i][1]].index);
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
        fprintf(gnu_files[3], "%lf %lf %d\n", center.x, center.y, i);
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
    fprintf(gnu_files[0], "plot '%s' using 1:2 with lines ls 1 title \"Tessellations\",", lines);
    fprintf(gnu_files[0], "'%s' using 1:2 with points pt 7 notitle,", datapoints);
    fprintf(gnu_files[0], "'' using 1:2:3 with labels point pt 7 offset char -1,-1 notitle,");
    fprintf(gnu_files[0], "'%s' using 1:2:3 with labels point pt 3 offset char -1,-1 notitle, ", extrapoints);
    fprintf(gnu_files[0], "'%s' using 1:2:3 with labels point pt 2 offset char -1,-1 notitle, ", centerpoint);
    fprintf(gnu_files[0], "'%s' using 1:2 with lines ls 2 title \"Calculated Path\", ", path);
    fprintf(gnu_files[0], "'%s' using 1:2 with lines ls 3 title \"Shortest Path\"\n", shortest);
    fclose(output);
    fclose(data);
    fclose(gnu_files[0]);
    fclose(gnu_files[1]);
    fclose(gnu_files[2]);
    fclose(gnu_files[3]);
    fclose(gnu_files[4]);
    fclose(gnu_files[5]);
    char plot[1024];
    sprintf(plot, "gnuplot -persistent %s", commands);
    system(plot);
    fclose(gnu_files[7]);
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
    /* calculate average point */
    for(i = 0; i < size; i++) {
        sum_x += points[i].x;
        sum_y += points[i].y;
    }
    center.x = sum_x / size;
    center.y = sum_y / size;
    /* plot center point */
    //fprintf(gnu_files[4], "%lf %lf %s\n", center.x, center.y, "C");
    fprintf(gnu_files[4], "%lf %lf\n", center.x, center.y);
    //printf("begin = %d\n", begin.index);
    //printf("\n");
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
        visited_count++;
    }
    return;
}

/* calculates the connections for un-joined vertices */
void join_vertex(vector<int *> *segments, struct point_t *points, struct point_t begin, int n, int size)
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
    }
    return;
}

/* calculates the connections for un-joined segments */
void join_segment(vector<int *> *segments, struct point_t *points, struct point_t begin, struct point_t end, int n, int m, int size)
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
    //printf("... joining: <%d,%d>\n", points[n].index, points[k].index);
    /* pushes new segments */
    pushed_segment = new int [2];
    pushed_segment[0] = tmp_segments[m][0];
    pushed_segment[1] = tmp_segments[m][1];
    segments->push_back(pushed_segment);
    return;
}

vector<struct polygon_t> construct_polygons(vector<int *> segments, struct point_t *points, int size)
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
vector<vector<int> > tessellate(vector<vector<int> > tessellations, vector<int *> segments, struct point_t *points, int size, char init, char add, int branch)
{
    vector<int *> edges; //pointer contains the index followed by position
    vector<int> *shape = NULL;
    struct vector_t X; //X-axis vector
    X.name = new char [2];
    X.name[0] = 'X';
    X.name[1] = '\0';
    struct vector_t Y; //Y-axis vector
    Y.name = new char [2];
    Y.name[0] = 'Y';
    Y.name[1] = '\0';
    struct point_t start;
    int i = 0;

    for(i = 0; i < size; i++) {
        /* find the initial cluster of edges */
        edges = edge_search(segments, points[i].index, points, size);
        /* skip if the requested branch is out of bounds */
        if(edges.size() < branch) {
            break;
        }
        start = points[i];
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
vector<struct polygon_t> tessellate_cross(vector<int *> segments, int i, int j, struct point_t *points, int size)
{
    vector<struct polygon_t> tessellations;
    struct polygon_t polygon;
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
    start = points[segments_i[i][0]];
    prev = points[segments_i[i][1]];
    Y.point[0].x = start.x;
    Y.point[0].y = start.y;
    Y.point[0].index = start.index;
    Y.i = (prev.x - start.x) / distance_p(prev, start);
    Y.j = (prev.y - start.y) / distance_p(prev, start);
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
    Y.point[0].x = start.x;
    Y.point[0].y = start.y;
    Y.point[0].index = start.index;
    Y.i = (prev.x - start.x) / distance_p(prev, start);
    Y.j = (prev.y - start.y) / distance_p(prev, start);
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
    Y.point[0].x = start.x;
    Y.point[0].y = start.y;
    Y.point[0].index = start.index;
    Y.i = (prev.x - start.x) / distance_p(prev, start);
    Y.j = (prev.y - start.y) / distance_p(prev, start);
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
    Y.point[0].x = start.x;
    Y.point[0].y = start.y;
    Y.point[0].index = start.index;
    Y.i = (prev.x - start.x) / distance_p(prev, start);
    Y.j = (prev.y - start.y) / distance_p(prev, start);
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
vector<int> *find_shape(vector<int *> segments, struct point_t *points, struct point_t start, int size, char init, char add, struct vector_t X, struct vector_t Y)
{
    struct point_t prev = start;
    struct point_t root = start;
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
vector<int> init_path(vector<int> path, vector<int *> edges, struct point_t *points, struct vector_t X, struct vector_t Y, char type, int branch)
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
    int count = -1;
    /* calculate and set flags for each edge */
    for(i = 0; i < edges.size(); i++) {
        /* set edge vector */
        E.point[0].x = points[edges[i][0]].x;
        E.point[0].y = points[edges[i][0]].y;
        E.point[0].index = points[edges[i][0]].index;
        E.point[1].x = points[edges[i][1]].x;
        E.point[1].y = points[edges[i][1]].y;
        E.point[1].index = points[edges[i][1]].index;
        E.i = (points[edges[i][1]].x - points[edges[i][0]].x);
        E.j = (points[edges[i][1]].y - points[edges[i][0]].y);
        E.length = length_v(E);
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
vector<struct polygon_t> optimize_polygons(vector<struct polygon_t> polygons, vector<int *> *segments, struct point_t *points, int size)
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
void remove_crosses(vector<int *> *segments, struct point_t *points, int size)
{
    struct polygon_t tmp_polygon;
    struct point_t p1;
    struct point_t p2;
    struct point_t p3;
    struct point_t p4;
    struct point_t tmp;
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

void finalize_segments(vector<int *> *segments, struct point_t *points, int size)
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
                if(intersection(points[i], points[j], points[(*segments)[k][0]], points[(*segments)[k][1]])) {
                    push = 0;
                    break;
                }
            }
            /* goes through all temporary segments */
            for(k = 0; k < tmp_segments.size(); k++) {
                if(intersection(points[i], points[j], points[tmp_segments[k][0]], points[tmp_segments[k][1]])) {
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
                if(intersection(points[i], points[j], points[(*segments)[k][0]], points[(*segments)[k][1]])) {
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

/* returns 1 if the points intersect, 0 if they don't */
int intersection(struct point_t p1, struct point_t p2, struct point_t p3, struct point_t p4)
{
    struct point_t tmp;
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
        return 0;
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
        return 0;
    }
    /* p3---------p4  p1---------p2 */
    if(l >= r) {
        return 0;
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
            return 1;
        }
        else if(p1.x < p4.x) {
            return 1;
        }
    }
    /* don't push the segment if there is an intersection */
    else if(x <= r && x >= l) {
        return 1;
    }
    return 0;
}

struct polygon_t find_shortest_path(vector<struct polygon_t> polygons, struct point_t *points, int size)
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
int accept_polygon(struct polygon_t polygon, vector<int *> segments, struct point_t *points) {
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

void visit_polygon(int *visited, struct polygon_t polygon, struct point_t *points)
{
    for(int i = 0; i < polygon.shape.size(); i++) {
        visited[points[polygon.shape[i]].index] = 1;
    }
}

struct polygon_t add_polygons(struct polygon_t A, struct polygon_t B, struct point_t *points)
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
struct polygon_t sub_polygons(struct polygon_t A, struct polygon_t B, struct point_t *points)
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
