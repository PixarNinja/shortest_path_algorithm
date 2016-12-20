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
#include <iterator>

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

/* global variables */
int permutations = 1;

void construct_segments(struct point_t *points, struct point_t begin, int n, int size, FILE *gnu_files[NUM_FILES], int *mapped, int **recorded, vector<int *> *segments);
vector<vector<int> > construct_polygons(vector<int *> segments, struct point_t *points, int size, FILE *gnu_files[NUM_FILES]);
vector<int *> edge_search(vector<int *> segments, int vertex);
vector<int> start_polygon(vector<int> path, vector<int *> edges, struct point_t *points, struct vector_t X, struct vector_t Y);
int shape_search(vector<int> shape, int vertex);
int polygons_search(vector<vector<int> > polygons, int vertex);
int segment_match(vector<int *> segments, int beginning, int end);
int duplicate_search(deque<int> queue, int vertex, int dups[2]);
deque<int> *separate_shape(deque<int> *queue, int size, int n, int m);
deque<int> *merge_queue(deque<int> *queue, int size);
double calculate_curvature(struct vector_t T1, struct vector_t T2, double tao);
double angle_t(double tao);
double tao_distance(struct vector_t V, double curvature, double theta);
double angle_v(struct vector_t V1, struct vector_t V2);
double distance_p(struct point_t start, struct point_t end);
double distance_v(struct vector_t V1, struct vector_t V2);
double length_v(struct vector_t V);
double dot_product(struct vector_t V1, struct vector_t V2);
void print_v(struct vector_t V);
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
            construct_segments(points, points[i], i, size, gnu_files, mapped, recorded, segments);
            keep_going = 0;
        }
        else {
            i = size;
        }
    }
    polygons = construct_polygons(*segments, points, size, gnu_files);
    /*
    printf("\n\nMAPPED ARRAY:\n");
    for(i = 0; i < size; i++) {
        printf("%d: %d\n", i, mapped[i]);
    }*/
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
vector<vector<int> > construct_polygons(vector<int *> segments, struct point_t *points, int size, FILE *gnu_files[NUM_FILES])
{
    struct vector_t X; //X-axis vector
    X.name = new char [2];
    X.name[0] = 'X';
    X.name[1] = '\0';
    struct vector_t Y; //Y-axis vector
    Y.name = new char [2];
    Y.name[0] = 'Y';
    Y.name[1] = '\0';
    struct vector_t S; //segment vector
    S.name = new char [2];
    S.name[0] = 'S';
    S.name[1] = '\0';
    struct point_t *search = new struct point_t [size];
    struct point_t prev;
    struct point_t curr;
    struct point_t start;
    struct point_t center;
    double sum_x = 0;
    double sum_y = 0;
    vector<vector<int> > polygons;
    vector<int *> edges; //contains the index followed by position
    vector<int> path; //contains the path of nodes to add as a polygon
    int i = 0;
    int j = 0;

    start = points[segments[0][0]];
    curr = start;
    /* calculate average point */
    for(i = 0; i < size; i++) {
        sum_x += points[i].x;
        sum_y += points[i].y;
    }
    center.x = sum_x / size;
    center.y = sum_y / size;
    /* initialize axis vectors */
    Y.point[0].x = start.x;
    Y.point[0].y = start.y;
    Y.point[0].index = start.index;
    Y.i = (start.x - center.x) / distance_p(center, start);
    Y.j = (start.y - center.y) / distance_p(center, start);
    Y.point[1].x = start.x + Y.i;
    Y.point[1].y = start.y + Y.j;
    Y.point[1].index = -1;
    Y.length = length_v(Y);
    print_v(Y);
    X.point[0].x = start.x;
    X.point[0].y = start.y;
    X.point[0].index = start.index;
    X.i = Y.j;
    X.j = -Y.i;
    X.point[1].x = start.x + X.i;
    X.point[1].y = start.y + X.j;
    X.point[1].index = -1;
    X.length = length_v(X);
    print_v(X);
    /* find the initial cluster of edges */
    edges = edge_search(segments, start.index);
    /* prints the initial edges */
    for(i = 0; i < edges.size(); i++) {
        printf("%d ", edges[i][1]);
    }
    /* find the initial direction */
    path = start_polygon(path, edges, points, X, Y);
    /* loop through path, always adding to the left */
    /* keep track of points to go back to later */
    /* loop until all points have added shapes */
    return polygons;
}

/* searches through a vector of segments for all matching end or
 * beginning segments, storing their indices in vector edges */
vector<int *> edge_search(vector<int *> segments, int vertex)
{
    vector<int *> edges;
    int *tmp;

    /* check if the segment is found */
    for(int i = 0; i < segments.size(); i++) {
        if(segments[i][0] == vertex) {
            tmp = new int [2];
            tmp[0] = vertex;
            tmp[1] = segments[i][1];
            edges.push_back(tmp);
        }
        else if(segments[i][1] == vertex) {
            tmp = new int [2];
            tmp[0] = vertex;
            tmp[1] = segments[i][0];
            edges.push_back(tmp);
        }
    }
    return edges;
}

/* finds the "right-most" path from the vertex */
vector<int> start_polygon(vector<int> path, vector<int *> edges, struct point_t *points, struct vector_t X, struct vector_t Y)
{
    int **X_quads = new int * [edges.size()];
    int **Y_quads = new int * [edges.size()];
    int i = 0;
    int j = 0;

    /* initialize the quadrant flags for each edge */
    for(i = 0; i < edges.size(); i++) {
        X_quads[i] = new int [4];
        Y_quads[i] = new int [4];
        for(j = 0; j < 4; j++) {
            X_quads[i][j] = -1;
            Y_quads[i][j] = -1;
        }
    }
    return path;
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

/* merges queue entries if they match */
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
