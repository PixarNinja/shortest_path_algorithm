#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <ctype.h>
#include <unistd.h>

#include <new>

#define NUM_FILES 4

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

double shortest_path(struct point_t *points, int size, FILE *gnu_files[NUM_FILES]);
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
    int size = 0;
    int permutations = 0;
    int i = 0;
    int j = 0;
    if(argc == 1) {
        printf("\n\nPlease enter the path of the .dat file to read from. Exiting Program. Good Day.\n\n");
    }
    gnu_files[0] = fopen ("./gnu_files/commands.tmp", "w+");
    gnu_files[1] = fopen("./gnu_files/points.tmp", "w+");
    gnu_files[2] = fopen("./gnu_files/lines.tmp", "w+");
    gnu_files[3] = fopen("./gnu_files/tmp.tmp", "w+");
    printf("%s\n", argv[1]);
    data = fopen(argv[1], "r");
    while(fgets(buf, 1024, data)) {
        size++;
    }
    fclose(data);
    points = new struct point_t [size];
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
    permutations += shortest_path(points, size, gnu_files);
    printf("\n");
    printf("Total Permutations: %d\n", permutations);
    printf("\n");
    /* plot */
    fclose(gnu_files[0]);
    fclose(gnu_files[1]);
    fclose(gnu_files[2]);
    system("gnuplot -persistent ./gnu_files/commands.tmp");
    fclose(gnu_files[3]);
    return 0;
}

/* calculates the shortest path */
double shortest_path(struct point_t *points, int size, FILE *gnu_files[NUM_FILES])
{
    struct vector_t V;
    struct vector_t T1;
    struct vector_t T2;
    struct vector_t initial;
    struct vector_t check;
    struct point_t *curr = new struct point_t [size];
    struct point_t *unvisited_points = new struct point_t [size];
    struct point_t *search = new struct point_t [size];
    struct point_t best;
    struct point_t begin;
    struct point_t start;
    struct point_t end;
    struct point_t prev;
    struct point_t center;
    double current_distance = DBL_MAX;
    double sum_x = 0.0;
    double sum_y = 0.0;
    double theta = DBL_MAX;
    double tmp = DBL_MAX;
    int **segments = new int * [size + 1];
    int ** recorded = new int * [size];
    int **loop = new int * [size];
    int *visited = new int [size];
    int *contoured_points = new int [size];
    int count = 0;
    int permutations = 0;
    int total_size = size;
    int index;
    int remaining_points = 0;
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;
    int m = 0;
    int n = 0;
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
        segments[i] = new int [2];
        segments[i][0] = INT_MAX;
        segments[i][1] = INT_MAX;
        recorded[i] = new int [size];
        for(j = 0; j < size; j++) {
            recorded[i][j] = 0;
        }
        loop[i] = new int [size];
    }
    segments[i] = new int [2];
    segments[i][0] = INT_MAX;
    segments[i][1] = INT_MAX;
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
    for(i = 0; i < size; i++) {
        if(distance_p(points[i], center) < tmp) {
            printf("data.[%d] < ", i);
            begin.x = points[i].x;
            begin.y = points[i].y;
            begin.index = points[i].index;
            tmp = distance_p(points[i], center);
        }
    }
    printf("\n");
    start.x = begin.x;
    start.y = begin.y;
    start.index = begin.index;
    end.x = begin.x;
    end.y = begin.y;
    end.index = begin.index;
    prev.x = begin.x;
    prev.y = begin.y;
    prev.index = begin.index;
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
    while(permutations <= total_size) {
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
            /* for debugging tao-distance function */
            print(V, T1, T2, curr[i].curvature, curr[i].theta, curr[i].tao);
            
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
            }
        }
        /* record segment */
        segments[k][0] = start.index;
        segments[k][1] = best.index;
        k++;
        /* record loops */
        loop[l][m++] = best.index;
        /* if the best point has been visited before */
        if(visited[best.index] == 1) {
            loop[l][m] = loop[l][1];
            /* prints the loop for debugging */
            printf("Loop %d: ", l);
            for(j = 0; j < m; j++) {
                printf("%d->", loop[l][j]);
            }
            printf("%d\n", loop[l][m]);
            /* checks if the loop has ended */
            if(begin.index == points[n].index) {
                permutations++;
                continue;
            }
            if(loop[l][0] == size) {
                /* don't reinitialize the initial vector */
            }
            else {
                /* initialize initial vector */
                for(j = 0; j < size; j++) {
                    if(loop[l][0] == search[j].index)
                        break;
                }
                initial.point[0].x = search[j].x;
                initial.point[0].y = search[j].y;
                initial.point[0].index = search[j].index;
                initial.i = (search[j].x - best.x) / distance_p(best, search[j]);
                initial.j = (search[j].y - best.x) / distance_p(best, search[j]);
                initial.point[1].x = search[j].x + initial.i;
                initial.point[1].y = search[j].y + initial.j;
                initial.point[1].index = size;
                initial.length = length_v(initial);
            }
            /* initializing vector V */
            V.point[0].x = search[j].x;
            V.point[0].y = search[j].y;
            V.point[0].index = search[j].index;
            V.i = 0;
            V.j = 0;
            V.length = 0;
            /* -- initializing vector T1 */
            T1.point[0].x = search[j].x;
            T1.point[0].y = search[j].y;
            T1.point[0].index = search[j].index;
            /* checks for the point with the smallest deviation in angle from
               the position vector of the starting point (selects the "second"
               point) */
            theta = DBL_MAX;
            for(j = 0; j < size; j++) {
                /* only check unvisited points */
                if(visited[search[j].index] == 0) {
                    check.point[0].x = start.x;
                    check.point[0].y = start.y;
                    check.point[0].index = start.index;
                    check.point[1].x = search[j].x;
                    check.point[1].y = search[j].y;
                    check.point[1].index = search[j].index;
                    check.i = check.point[1].x - check.point[0].x;
                    check.j = check.point[1].y - check.point[0].y;
                    check.length = length_v(check);
                    tmp = angle_v(T2, check);
                    if(tmp < theta) {
                        theta = tmp;
                        index = j;
                    }
                }
            }
            /* points T1 in the direction of the "second point" */
            T1.i = (search[index].x - start.x) / distance_p(search[index], start);
            T1.j = (search[index].y - start.y) / distance_p(search[index], start);
            T1.point[1].x = start.x + T1.i;
            T1.point[1].y = start.y + T1.j;
            T1.point[1].index = INT_MAX;
            T1.length = 1;
            /* initializing vector T2 */
            T2.point[0].x = start.x;
            T2.point[0].y = start.y;
            T2.point[0].index = start.index;
            /* stores index of the starting point */
            index = start.index;
            count = 0;
            permutations++;
            m = 0;
            loop[l++][m++] = size;
            continue;
        }
        /* if the loop is compete TODO: finish implementation */
        if(0) {
            /* finds starting point of the next curve */
            remaining_points = 0;
            for(j = 0; j < size; j++) {
                if(visited[j] == 0) {
                    unvisited_points[remaining_points++] = search[j];
                }
                else {
                    contoured_points[j] = 1;
                }
            }
            current_distance = DBL_MAX;
            if(remaining_points == 0) {
                count = size;
                permutations = total_size + 1;
                continue;
            }
            for(j = 0; j < remaining_points; j++) {
                if(distance_p(best, unvisited_points[j]) < current_distance) {
                    current_distance = distance_p(start, unvisited_points[j]);
                    start = unvisited_points[j];
                }
            }
            /* sets new ending point */
            end = start;
            /* calculate new average point */
            sum_x = 0;
            sum_y = 0;
            center.x = sum_x / remaining_points;
            center.y = sum_y / remaining_points;
            visited[start.index] = 1;
            /* initialize initial vector */
            initial.point[0].x = start.x;
            initial.point[0].y = start.y;
            initial.point[0].index = start.index;
            initial.i = (start.x - center.x) / distance_p(start, center);
            initial.j = (start.y - center.y) / distance_p(start, center);
            initial.point[1].x = start.x + initial.i;
            initial.point[1].y = start.y + initial.j;
            initial.point[1].index = INT_MAX;
            initial.length = length_v(initial);
            /* initializing point 1 data for structure k
            -- initializing vector V */
            V.point[0].x = start.x;
            V.point[0].y = start.y;
            V.point[0].index = start.index;
            V.i = 0;
            V.j = 0;
            V.length = 0;
            /* -- initializing vector T1 */
            T1.point[0].x = start.x;
            T1.point[0].y = start.y;
            T1.point[0].index = start.index;
            /* checks for the point with the smallest deviation in angle from
               the position vector of the starting point (selects the "second"
               point) */
            theta = DBL_MAX;
            for(j = 0; j < size; j++) {
                /* skip the starting point and the previous point */
                if((search[j].index != start.index) || (search[j].index != prev.index)) {
                    /* only check unvisited points */
                    if(visited[search[j].index] == 0) {
                        check.point[0].x = start.x;
                        check.point[0].y = start.y;
                        check.point[0].index = start.index;
                        check.point[1].x = search[j].x;
                        check.point[1].y = search[j].y;
                        check.point[1].index = search[j].index;
                        check.i = check.point[1].x - check.point[0].x;
                        check.j = check.point[1].y - check.point[0].y;
                        check.length = length_v(check);
                        tmp = angle_v(initial, check);
                        if(tmp < theta) {
                            theta = tmp;
                            index = j;
                        }
                    }
                }
            }
            /* points T1 in the direction of the "second point" */
            T1.i = (search[index].x - start.x) / distance_p(search[index], start);
            T1.j = (search[index].y - start.y) / distance_p(search[index], start);
            T1.point[1].x = start.x + T1.i;
            T1.point[1].y = start.y + T1.j;
            T1.point[1].index = INT_MAX;
            T1.length = 1;
            /* -- initializing vector T2 */
            T2.point[0].x = start.x;
            T2.point[0].y = start.y;
            T2.point[0].index = start.index;
            /* stores index of the starting point */
            index = start.index;
            count = 0;
            permutations++;
            continue;
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
        start = best;
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
    /* calculates the all segments for each contour */
    for(j = 0; j < k; j++) {
        /* skips over already-recorded segments */
        if(recorded[segments[j][0]][segments[j][1]] == 1) {
            continue;   
        }
        fprintf(gnu_files[2], "%lf %lf %d\n", points[segments[j][0]].x, points[segments[j][0]].y, points[segments[j][0]].index);
        fprintf(gnu_files[2], "%lf %lf %d\n", points[segments[j][1]].x, points[segments[j][1]].y, points[segments[j][1]].index);
        fprintf(gnu_files[2], "\n");
        recorded[segments[j][0]][segments[j][1]] = 1;
        /* book-keeping */
        printf("%d = (%d, %d): <%d,%d>\n", j, segments[j][0], segments[j][1], points[segments[j][0]].index, points[segments[j][1]].index);
    }
    /* final point
    fprintf(gnu_files[2], "%lf %lf %d\n", end.x, end.y, end.index);*/
    /* plot */
    fprintf(gnu_files[0], "plot './gnu_files/lines.tmp' using 1:2 with lines ls 1 title \"shortest path\",");
    fprintf(gnu_files[0], "'./gnu_files/points.tmp' using 1:2 with points pt 7 notitle,");
    fprintf(gnu_files[0], "'' using 1:2:3 with labels point pt 7 offset char -1,-1 notitle\n");
    return permutations;
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
