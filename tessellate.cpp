#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <ctype.h>
#include <unistd.h>

#include <new>
#include <vector>
#include <deque>
#include <algorithm>

#define NUM_FILES 4

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
    struct point_t point[2];
    double length;
    double i;
    double j;
};

/* global variables */
int permutations = 1;

void construct_segments(struct point_t *points, struct point_t begin, int n, int size, FILE *gnu_files[NUM_FILES], int *mapped, int **recorded, vector<int *> *segments);
vector<vector<int> > construct_polygons(vector<int *> segments, int size);
int shape_search(vector<int> shape, int vertex);
int polygons_search(vector<vector<int> > polygons, int vertex);
int segment_search(vector<int *> segments, int vertex);
int duplicate_search(deque<int> queue, int vertex, int dups[2]);
deque<int> *separate_shape(deque<int> *queue, int size, int n, int m);
deque<int> *merge_queue(deque<int> *queue, int size);
double calculate_curvature(struct vector_t T1, struct vector_t T2, double tao);
double calculate_theta(double tao);
double tao_distance(struct vector_t V, double curvature, double theta);
double angle_v(struct vector_t V1, struct vector_t V2);
double distance_p(struct point_t start, struct point_t end);
double distance_v(struct vector_t V1, struct vector_t V2);
double length_v(struct vector_t V);
double dot_product(struct vector_t V1, struct vector_t V2);
void print(struct vector_t V, struct vector_t T1, struct vector_t T2, double curvature, double theta, double tao);
void memory_error(void);

int main(int argc, char *argv[])
{
    FILE *data;
    FILE *gnu_files[NUM_FILES];
    struct point_t *points;
    char buf[1024];
    double range = 0.0;
    vector<vector<int> > polygons;
    vector<int *> *segments = new vector<int *> [1];
    int **recorded;
    int *mapped;
    int keep_going = 0;
    int size = 0;
    int i = 0;
    int j = 0;
    if(argc == 1) {
        printf("\n\nPlease enter the path of the .dat file to read from. Exiting Program. Good Day.\n\n");
    }
    gnu_files[0] = fopen ("./gnu_files/commands.polygons", "w+");
    gnu_files[1] = fopen("./gnu_files/points.polygons", "w+");
    gnu_files[2] = fopen("./gnu_files/lines.polygons", "w+");
    gnu_files[3] = fopen("./gnu_files/polygons.polygons", "w+");
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
    /* stores data for gnu_points */
    for(i = 0; i < size; i++) {
        fprintf(gnu_files[1], "%lf %lf %d\n", points[i].x, points[i].y, points[i].index);
    }
    /* plot setup */
    fprintf(gnu_files[0], "set xrange [%lf:%lf]\n", -(range + 1), range + 1);
    fprintf(gnu_files[0], "set yrange [%lf:%lf]\n", -(range + 1), range + 1);
    fprintf(gnu_files[0], "set size ratio 1\n");
    fprintf(gnu_files[0], "set grid\n");
    fprintf(gnu_files[0], "set title \"Contour Construction Algorithm\"\n");
    fprintf(gnu_files[0], "set style line 1 lc rgb \"black\" lw 1\n");
    /* runs tao-distance algorithm on data */
    for(i = 0; i < size; i++) {
        for(j = 0; j < size; j++) {
            if(mapped[j] == 0) {
                keep_going = 1;
            }
        }
        if(keep_going == 1) {
            construct_segments(points, points[i], i, size, gnu_files, mapped, recorded, segments);
            keep_going = 0;
        }
        else {
            i = size;
        }
    }
    polygons = construct_polygons(*segments, size);
    printf("\n\nMAPPED ARRAY:\n");
    for(i = 0; i < size; i++) {
        printf("%d: %d\n", i, mapped[i]);
    }
    /* print segment information */
    printf("\nSEGMENTS:\n");
    for(i = 0; i < segments->size(); i++) {
        printf("%d: <%d,%d>\n", i, (*segments)[i][0], (*segments)[i][1]);
    }
    /* print polygon information */
    printf("\nPOLYGONS:\n");
    for(i = 0; i < polygons.size(); i++) {
        printf("%d: ", i);
        for(j = 0; j < polygons[i].size(); j++) {
            printf("%d ", polygons[i][j]);
        }
        printf("\n");
    }
    /* plot */
    fprintf(gnu_files[0], "plot './gnu_files/lines.polygons' using 1:2 with lines ls 1 title \"shortest path\",");
    fprintf(gnu_files[0], "'./gnu_files/points.polygons' using 1:2 with points pt 7 notitle,");
    fprintf(gnu_files[0], "'' using 1:2:3 with labels point pt 7 offset char -1,-1 notitle\n");
    printf("\n");
    printf("Total Permutations: %d\n", permutations);
    printf("\n");
    fclose(gnu_files[0]);
    fclose(gnu_files[1]);
    fclose(gnu_files[2]);
    system("gnuplot -persistent ./gnu_files/commands.polygons");
    fclose(gnu_files[3]);
    return 0;
}

/* calculates the shortest path */
void construct_segments(struct point_t *points, struct point_t begin, int n, int size, FILE *gnu_files[NUM_FILES], int *mapped, int **recorded, vector<int *> *segments)
{
    struct vector_t V;
    struct vector_t T1;
    struct vector_t T2;
    struct point_t *curr = new struct point_t [size];
    struct point_t *search = new struct point_t [size];
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
    int segment_count = 0;
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
        search[i].x = points[i].x;
        search[i].y = points[i].y;
        search[i].tao_distance = points[i].tao_distance;
        search[i].index = points[i].index;
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
        /* store start index in visted-array */
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
            if((search[i].index == best.index) || (search[i].index == prev.index)) {
                curr[i].tao_distance = DBL_MAX;
                i++;
                count++;
                continue;
            }
            /* initializing vector V */
            V.point[1].x = search[i].x;
            V.point[1].y = search[i].y;
            V.point[1].index = search[i].index;
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
            curr[i].theta = calculate_theta(curr[i].tao);
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
                /* check if the loop has finished elsewhere */
                if(loop[m] != loop[j]) {
                    for(j = m; j > 0; j--) {
                        if(loop[j] != loop[0]) {
                            m--;
                        }
                        else {
                            break;
                        }
                    }
                }
                for(j = 0; j < m; j++) {
                    printf("%d->", search[loop[j]].index);
                    mapped[loop[j]] = 1;
                }
                printf("%d\n", search[loop[m]].index);
                for(j = 0; j < m; j++) {
                    /* record segment */
                    tmp_segments[j][0] = loop[j];
                    tmp_segments[j][1] = loop[j + 1];
                }
                tmp_segments[m][0] = loop[m];
                tmp_segments[m][1] = loop[0];
                /* calculates the tmp_segments for each contour */
                for(j = 0; j < m; j++) {
                    /* skips over recorded tmp_segments */
                    if((recorded[tmp_segments[j][0]][tmp_segments[j][1]] == 1) || (recorded[tmp_segments[j][1]][tmp_segments[j][0]] == 1)) {
                        continue;   
                    }
                    /* pushes new segments */
                    pushed_segment = new int [2];
                    pushed_segment[0] = points[tmp_segments[j][0]].index;
                    pushed_segment[1] = points[tmp_segments[j][1]].index;
                    segments->push_back(pushed_segment);
                    segment_count++;
                    /* plot */
                    fprintf(gnu_files[2], "%lf %lf %d\n", points[tmp_segments[j][0]].x, points[tmp_segments[j][0]].y, points[tmp_segments[j][0]].index);
                    fprintf(gnu_files[2], "%lf %lf %d\n", points[tmp_segments[j][1]].x, points[tmp_segments[j][1]].y, points[tmp_segments[j][1]].index);
                    fprintf(gnu_files[2], "\n");
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

/* calculate polygons given all contoured segments */
vector<vector<int> > construct_polygons(vector<int *> segments, int size)
{
    vector<vector<int> > polygons;
    vector<int *> tmp_segments;
    vector<int *> free_segments; //list of segments that are available for adding
    deque<int> *queue = new deque<int> [size]; //keeps track of adding and splitting
    vector<int> tmp_shape_0;
    vector<int> tmp_shape_1;
    vector<int> tmp_shape_2;
    int *pushed_segment;
    int found_beginning[2];
    int found_end[2];
    int dups[2] = {INT_MAX};
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;
    int m = 0;
    int n = 0;
    int start = size;
    int stop = size;
    int split = 0;
    int add = 0;

    tmp_shape_0.push_back(segments[0][0]);
    tmp_shape_0.push_back(segments[0][1]);
    free_segments.push_back(segments[0]);
    while(segments[1][1] != segments[0][0]) {
        tmp_shape_0.push_back(segments[1][1]);
        free_segments.push_back(segments[1]);
        segments.erase(segments.begin() + 1);
    }
    free_segments.push_back(segments[1]);
    printf("\nFREE:\n");
    for(i = 0; i < free_segments.size(); i++) {
        printf("<%d,%d>\n", free_segments[i][0], free_segments[i][1]);
    }
    segments.erase(segments.begin(), segments.begin() + 2);
    //segments.erase(segments.begin() + 1);
    polygons.push_back(tmp_shape_0);
    /* reorder the segments closest to the first polygon */
    printf("\nNEW SEGMENTS:\n");
    for(i = 0; i < segments.size(); i++) {
        /* add the segment if there is a matching vertex */
        if(shape_search(polygons[0], segments[i][0]) || shape_search(polygons[0], segments[i][1])) {
            pushed_segment = new int [2];
            pushed_segment[0] = segments[i][0];
            pushed_segment[1] = segments[i][1];
            tmp_segments.push_back(pushed_segment);
            segments.erase(segments.begin() + i);
            i--;
            printf("<%d,%d>\n", pushed_segment[0], pushed_segment[1]);
        }
    }
    /* reorder the segments closest to the reordered segments */
    for(i = 0; i < segments.size(); i++) {
        /* add the segment if there is a matching vertex */
        if((segment_search(tmp_segments, segments[i][0]) != INT_MAX) || (segment_search(tmp_segments, segments[i][1]) != INT_MAX)) {
            pushed_segment = new int [2];
            pushed_segment[0] = segments[i][0];
            pushed_segment[1] = segments[i][1];
            tmp_segments.push_back(pushed_segment);
            segments.erase(segments.begin() + i);
            i--;
            printf("<%d,%d>\n", pushed_segment[0], pushed_segment[1]);
        }
    }
    printf("\n");
    printf("\nPOLYGON [0]:\n");
    for(j = 0; j < polygons[0].size(); j++) {
        printf("%d ", polygons[0][j]);
    }
    printf("\n");
    printf("\n");
    segments = tmp_segments;
    /* loop through all segments */
    for(i = 0; i < segments.size(); i++) {
        found_beginning[0] = INT_MAX;
        found_end[0] = INT_MAX;
        split = 0;
        /* check for a split */
        for(k = 0; k < polygons.size(); k++) {
            if(shape_search(polygons[k], segments[i][0]) && shape_search(polygons[k], segments[i][1])) {
                printf("SPLIT %d with <%d,%d>\n", k, segments[i][0], segments[i][1]);
                tmp_shape_1.clear();
                tmp_shape_2.clear();
                l = distance(polygons[k].begin(), find(polygons[k].begin(), polygons[k].end(), segments[i][0]));
                start = segments[i][0];
                stop = segments[i][1];
                tmp_shape_1.push_back(polygons[k][l]);
                /* leftward loop */
                while(start != stop) {
                    l--;
                    if(l < 0) {
                        l = polygons[k].size() - 1;
                    }
                    tmp_shape_1.push_back(polygons[k][l]);
                    start = polygons[k][l];
                }
                l = distance(polygons[k].begin(), find(polygons[k].begin(), polygons[k].end(), segments[i][0]));
                start = segments[i][0];
                stop = segments[i][1];
                tmp_shape_2.push_back(polygons[k][l]);
                /* rightward loop */
                while(start != stop) {
                    l++;
                    if(l > polygons[k].size() - 1) {
                        l = 0;
                    }
                    tmp_shape_2.push_back(polygons[k][l]);
                    start = polygons[k][l];
                }
                polygons.push_back(tmp_shape_1);
                polygons.push_back(tmp_shape_2);
                polygons.erase(polygons.begin() + k);
                split = 1;
                k += 2;
            }
            else if(shape_search(polygons[k], segments[i][0])) {
                found_beginning[0] = k;
                found_beginning[1] = segments[i][0];
            }
            else if(shape_search(polygons[k], segments[i][1])) {
                found_end[0] = k;
                found_end[1] = segments[i][1];
            }
        }
        if((found_beginning[0] != INT_MAX) && (found_end[0] != INT_MAX)) {
            printf("SPLIT %d and %d with <%d,%d>\n", found_beginning[0], found_end[0], found_beginning[1], found_end[1]);
            split = 1;
        }
        if(split) {
            continue;
        }
        /* update the queue */
        for(j = 0; j < size; j++) {
            /* add segment if the slot is empty */
            if(queue[j].size() < 1) {
                queue[j].push_back(segments[i][0]);
                queue[j].push_back(segments[i][1]);
                printf("add <%d,%d>\n", segments[i][0], segments[i][1]);
                break;
            }
            /* add segments[i][1] to the front of the list */
            if(segments[i][0] == queue[j][0]) {
                queue[j].push_front(segments[i][1]);
                printf("add %d to front\n", segments[i][1]);
                break;
            }
            /* add segments[i][0] to the front of the list */
            if(segments[i][1] == queue[j][0]) {
                queue[j].push_front(segments[i][0]);
                printf("add %d to front\n", segments[i][0]);
                break;
            }
            /* add segments[i][1] to the back of the list */
            if(segments[i][0] == queue[j][queue[j].size() - 1]) {
                queue[j].push_back(segments[i][1]);
                printf("add %d to back\n", segments[i][1]);
                break;
            }
            /* add segments[i][0] to the back of the list */
            if(segments[i][1] == queue[j][queue[j].size() - 1]) {
                queue[j].push_back(segments[i][0]);
                printf("add %d to back\n", segments[i][0]);
                break;
            }
        }
        /* merges queue elements if ends have a matching node */
        queue = merge_queue(queue, size);
        /* check for an addition */
        for(j = 0; j < size; j++) {
            if(queue[j].size() <= 1) {
                continue;
            }
            for(k = 0; k < polygons.size(); k++) {
                /* check if the addition is complete */
                if(shape_search(polygons[k], queue[j][0]) && shape_search(polygons[k], queue[j][queue[j].size() - 1])) {
                    printf("ADD to %d: <%d,%d>\n", k, queue[j][0], queue[j][queue[j].size() - 1]);
                    add = 1;
                    /* add the new shape to free_segments */
                    for(l = 0; l < queue[j].size() - 1; l++) {
                        pushed_segment = new int [2];
                        pushed_segment[0] = queue[j][l];
                        pushed_segment[1] = queue[j][l + 1];
                        free_segments.push_back(pushed_segment);
                    }
                    l = segment_search(free_segments, queue[j][queue[j].size() - 1]);
                    free_segments.erase(free_segments.begin() + l);
                    printf("\nFREE:\n");
                    for(l = 0; l < free_segments.size(); l++) {
                        printf("<%d,%d>\n", free_segments[l][0], free_segments[l][1]);
                    }
                    /* add the new shape */
                    tmp_shape_0.clear();
                    for(l = 0; l < queue[j].size(); l++) {
                        tmp_shape_0.push_back(queue[j][l]);
                    }
                    polygons.push_back(tmp_shape_0);
                    queue[j].erase(queue[j].begin(), queue[j].end());
                    for(; j < size - 1; j++) {
                        queue[j] = queue[j + 1];
                    }
                    queue[size - 1].clear();
                    break;
                }
            }
            /* check if the addition is a loop */
            for(k = 0; k < queue[j].size(); k++) {
                if(duplicate_search(queue[j], queue[j][k], dups)) {
                    printf("ADD loop around %d\n", queue[j][k], queue[j][k]);
                    add = 1;
                    /* add the new segments to free_segments */
                    n = dups[0];
                    m = dups[1];
                    pushed_segment = new int [2];
                    pushed_segment[0] = queue[j][n];
                    pushed_segment[1] = queue[j][n + 1];
                    free_segments.push_back(pushed_segment);
                    for(; n < m; n++) {
                        pushed_segment = new int [2];
                        pushed_segment[0] = queue[j][n];
                        pushed_segment[1] = queue[j][n + 1];
                        free_segments.push_back(pushed_segment);
                    }
                    l = segment_search(free_segments, queue[j][m]);
                    free_segments.erase(free_segments.begin() + l);
                    printf("\nFREE:\n");
                    for(l = 0; l < free_segments.size(); l++) {
                        printf("<%d,%d>\n", free_segments[l][0], free_segments[l][1]);
                    }
                    queue = separate_shape(queue, size, j, k);
                    tmp_shape_1.clear();
                    l = distance(queue[j].begin(), find(queue[j].begin(), queue[j].end(), queue[j][k]));
                    start = INT_MAX;
                    stop = queue[j][k];
                    /* adding to a temporary shape */
                    while(start != stop) {
                        l++;
                        if(l > queue[j].size() - 1) {
                            l = 0;
                        }
                        tmp_shape_1.push_back(queue[j][l]);
                        start = queue[j][l];
                        queue[j].erase(queue[j].begin() + l);
                        l--;
                    }
                    polygons.push_back(tmp_shape_1);
                }
            }
        }
        if(add) {
            continue;
        }
        /* check for an addition accross multiple shapes */
        for(j = 0; j < size; j++) {
            if(queue[j].size() <= 1) {
                continue;
            }
            if(((n = polygons_search(polygons, queue[j][0])) != INT_MAX) && ((m = polygons_search(polygons, queue[j][queue[j].size() - 1])) != INT_MAX)) {
                if(n == m) {
                    continue;
                }
                printf("ADD to %d and %d\n", n, m);
                /* use free_segments to find if a polygon exists between the points */
                n = queue[j][0];
                m = queue[j][queue[j].size() - 1];
                while(n != m) {
                    /* if there is a match in find_segments */
                    if((l = segment_search(free_segments, n)) != INT_MAX) {
                        n = l;
                    }
                    else {
                        printf("\nNO POLYGON\n");
                        break;
                    }
                }
            }
        }
    }
    printf("\nQUEUE:\n");
    for(i = 0; i < size; i++) {
        if(queue[i].size() >= 1) {
            printf("%d: ", i);
            for(j = 0; j < queue[i].size(); j++) {
                printf("%d ", queue[i][j]);
            }
            printf("\n");
        }
    }
    printf("\nFREE:\n");
    for(i = 0; i < free_segments.size(); i++) {
        printf("<%d,%d>\n", free_segments[i][0], free_segments[i][1]);
    }
    return polygons;
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
    return INT_MAX;
}

/* searches through a vector of segments for a matching vertex */
int segment_search(vector<int *> segments, int vertex)
{
    /* check if the segment is found */
    for(int i = 0; i < segments.size(); i++) {
        if((segments[i][0] == vertex) || (segments[i][1] == vertex)) {
            return i;
        }
    }
    return INT_MAX;
}

/* searches through a deque of vertices for a duplicate */
int duplicate_search(deque<int> queue, int vertex, int dups[2])
{
    int d1 = 0;
    int d2 = 0;
    /* check if the vertex has a duplicate */
    if(find(queue.begin(), queue.end(), vertex) != queue.end()) {
        d1 = distance(queue.begin(), find(queue.begin(), queue.end(), vertex));
        queue.erase(queue.begin(), queue.begin() + d1 + 1);
        if(find(queue.begin() + d1, queue.end(), vertex) != queue.end()) {
            d2 = distance(queue.begin(), find(queue.begin(), queue.end(), vertex));
            dups[0] = d1;
            dups[1] = d1 + d2;
            return 1;
        }
    }
    return 0;
}

/* separates a shape given there is a found duplicate at index m */
deque<int> *separate_shape(deque<int> *queue, int size, int n, int m)
{
    int i = 0;
    int j = 0;
    for(i = 0; i < size; i++) {
        /* add segment to an empty slot */
        if(queue[i].size() < 1) {
            for(j = m + 1; queue[n][j] != queue[n][m]; j++) {
                queue[i].push_back(queue[n][j]);
            }
            queue[i].push_back(queue[n][m]);
            break;
        }
    }
    return queue;
}

deque<int> *merge_queue(deque<int> *queue, int size)
{
    int i = 0;
    int j = 0;
    int k = 0;
    for(i = 0; i < size; i++) {
        for(j = 0; j < size; j++) {
            if((i == j) || (queue[i].size() < 1) || (queue[j].size() < 1)) {
                continue;
            }
            /* check if queue[i][last] matches queue[j][0] */
            if((queue[i][queue[i].size() - 1]) == queue[j][0]) {
                printf("MERGE\n");
                for(k = 1; k < queue[j].size(); k++) {
                    queue[i].push_back(queue[j][k]);
                }
                queue[j].erase(queue[j].begin(), queue[j].end());
                for(k = j; k < size - 1; k++) {
                    queue[k] = queue[k + 1];
                }
                queue[size - 1].clear();
                break;
            }
        }
    }
    return queue;
}

/* calculates curvature given structure k */
double calculate_curvature(struct vector_t T1, struct vector_t T2, double tao)
{
    return (distance_v(T1, T2) / calculate_theta(tao));
}

/* calculates theta given structure k */
double calculate_theta(double tao)
{
    return (acos(tao) + (M_PI / 180));
    //return (acos(tao));
}

/* calculates distance given index and structure */
double tao_distance(struct vector_t V, double curvature, double theta)
{
    return (V.length + curvature + theta);
    //return (V.length + curvature + (cos(theta) * V.length));
}

/* calculates angle between two vectors */
double angle_v(struct vector_t V1, struct vector_t V2)
{
    return (acos(dot_product(V1, V2) / (V1.length * V2.length)));
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

/* prints structure k for debugging */
void print(struct vector_t V, struct vector_t T1, struct vector_t T2, double curvature, double theta, double tao)
{
    printf("V: %d[0](%lf, %lf), %d[1](%lf, %lf), <%lf, %lf>, |V| = %lf\n", V.point[0].index, V.point[0].x, V.point[0].y, V.point[1].index, V.point[1].x, V.point[1].y, V.i, V.j, V.length);
    printf("T1: point[0](%lf, %lf), point[1](%lf, %lf), <%lf, %lf>, |T1| = %lf\n", T1.point[0].x, T1.point[0].y, T1.point[1].x, T1.point[1].y, T1.i, T1.j, T1.length);
    printf("T2: point[0](%lf, %lf), point[1](%lf, %lf), <%lf, %lf>, |T2| = %lf\n", T2.point[0].x, T2.point[0].y, T2.point[1].x, T2.point[1].y, T2.i, T2.j, T2.length);
    printf("curvature: %lf; ", curvature);
    printf("angle = %lf; ", theta * 180 / M_PI);
    printf("tao = %lf; ", tao);
    printf("tao-distance = %lf\n\n", V.point[1].tao_distance);
}

/* prints to the terminal if there is an error assigning memory */
void memory_error(void)
{
    printf("\n\nError assigning memory. Exiting Program. Good Day.\n\n");
}
