#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "construct_contour.h"

#define FILE_NUM 6

struct function_t {
    struct point_t point1;
    struct point_t point2;
    double m;
};

struct triangle_t {
    char name[3];
    char *points;
    struct function_t f1;
    struct function_t f2;
    struct function_t f3;
    double x;
    double y;
};

int global_count = 0;

double path(struct point_t *point, int size, FILE *gnu_files[FILE_NUM]);
void construct_neighbors(int **neighbors, struct point_t *point, int size, FILE *gnu_files[FILE_NUM]);
struct triangle_t *construct_triangles(struct point_t **curve, int size, int divisions);

int main(void)
{
    FILE *gnu_files[FILE_NUM];
    gnu_files[0] = fopen ("./gnu_files/commands.tmp", "w+");
    gnu_files[1] = fopen("./gnu_files/points.tmp", "w+");
    gnu_files[2] = fopen("./gnu_files/lines.tmp", "w+");
    gnu_files[3] = fopen("./gnu_files/vectors_t.tmp", "w+");
    gnu_files[4] = fopen("./gnu_files/vectors_n.tmp", "w+");
    gnu_files[5] = fopen("./gnu_files/tmp.tmp", "w+");
    FILE *data = fopen("./datapoints/datapoints.dat", "r");
    struct point_t *point;
    char buf[1024];
    double distance = 0;
    int size = 0;
    int i = 0;
    while(fgets(buf, 1024, data)) {
        size++;
    }
    fclose(data);
    data = fopen("./datapoints/datapoints.dat", "r");
    point = malloc(sizeof(struct point_t) * size);
    while(fscanf(data, "%d: (%lf, %lf)", &point[i].index, &point[i].x, &point[i].y) > 0) {
        i++;
    }
    fclose(data);
    /* stores data for gnu_points */
    for(i = 0; i < size; i++) {
        fprintf(gnu_files[1], "%lf %lf %d\n", point[i].x, point[i].y, point[i].index);
    }
    /* plot setup */
    fprintf(gnu_files[0], "set xrange [-8:8]\n");
    fprintf(gnu_files[0], "set yrange [-8:8]\n");
    fprintf(gnu_files[0], "set size ratio 1\n");
    fprintf(gnu_files[0], "set grid\n");
    fprintf(gnu_files[0], "set title \"DATA POINTS\"\n");
    fprintf(gnu_files[0], "set style line 1 lc rgb \"black\" lw 1\n");
    fprintf(gnu_files[0], "set style line 2 lc rgb \"grey\"\n");
    fprintf(gnu_files[0], "set style line 3 lc rgb \"orange\"\n");
    fprintf(gnu_files[0], "set style line 4 lc rgb \"cyan\"\n");
    fprintf(gnu_files[0], "set style arrow 1 head filled size screen 0.025,20,35 ls 3\n");
    fprintf(gnu_files[0], "set style arrow 2 head filled size screen 0.025,20,35 ls 4\n");
    /* runs tao-distance algorithm on data */
    distance = path(point, size - 1, gnu_files);
    printf("\n");
    printf("Total Permutations: %d\n\n", global_count);
    printf("Shortest Path: %s\n", "derp");
    printf("Distance: %lf\n", distance);
    printf("\n");
    fclose(gnu_files[0]);
    fclose(gnu_files[1]);
    fclose(gnu_files[2]);
    fclose(gnu_files[3]);
    fclose(gnu_files[4]);
    system("gnuplot -persistent ./gnu_files/commands.tmp");
    return 0;
}

/* calculates the shortest path */
double path(struct point_t *point, int size, FILE *gnu_files[FILE_NUM])
{
    double total = size;
    int **neighbors = malloc(sizeof(int *) * (size + 1));
    int i = 0;
    for(i = 0; i <= size; i++) {
        /* each point can have 0, 1, or 2 neighbors */
        neighbors[i] = malloc(sizeof(int) * 2);
        neighbors[i][0] = size + 1;
        neighbors[i][1] = size + 1;
    }
    construct_neighbors(neighbors, point, size, gnu_files);
    /* printing for debug */
    printf("NEIGHBORS ARRAY:\n");
    for(i = 0; i <= size; i++) {
        printf("node %d: %d, %d\n", i, neighbors[i][0], neighbors[i][1]);
        if(neighbors[i][0] == (size + 1)) {
            fprintf(gnu_files[2], "%lf %lf\n", point[neighbors[i][1]].x, point[neighbors[i][1]].y);
            fprintf(gnu_files[2], "%lf %lf\n", point[i].x, point[i].y);
            fprintf(gnu_files[2], "%lf %lf\n", point[neighbors[i][1]].x, point[neighbors[i][1]].y);
        }
        else {
            fprintf(gnu_files[2], "%lf %lf\n", point[neighbors[i][0]].x, point[neighbors[i][0]].y);
            fprintf(gnu_files[2], "%lf %lf\n", point[i].x, point[i].y);
            fprintf(gnu_files[2], "%lf %lf\n", point[neighbors[i][1]].x, point[neighbors[i][1]].y);
        }
        fprintf(gnu_files[2], "\n");
    }
    return total;
}

/* creates a 2D array (size x 2) of neighbors from the point data */
void construct_neighbors(int **neighbors, struct point_t *point, int size, FILE *gnu_files[FILE_NUM])
{
    struct point_t *contour_point = malloc(sizeof(struct point_t) * (size + 1));
    struct point_t center;
    struct vector_t *position = malloc(sizeof(struct vector_t) * (size + 1));
    double *distance = malloc(sizeof(double) * (size + 1));
    double *circles = malloc(sizeof(double) * (size + 1));
    double range[2] = {0};
    double division_constant = 0;
    double shortest = DBL_MAX;
    double longest = 0;
    double tao = 0;
    double sum_x = 0;
    double sum_y = 0;
    int **curve = malloc(sizeof(int *) * (size + 1));
    int **tmp = malloc(sizeof(struct point_t) * (size + 1));
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;
    int m = 0;
    int division_number = 0;
    /* calculate average point */
    for(i = 0; i <= size; i++) {
        sum_x += point[i].x;
        sum_y += point[i].y;
    }
    center.x = sum_x / size;
    center.y = sum_y / size;
    fprintf(gnu_files[1], "%lf %lf %s\n", center.x, center.y, "C");
    /* calculate position vectors */
    for(i = 0; i <= size; i++) {
        position[i].point[0].x = center.x;
        position[i].point[0].y = center.y;
        position[i].point[0].index = size + 1;
        position[i].point[1].x = point[i].x;
        position[i].point[1].y = point[i].y;
        position[i].point[1].index = point[i].index;
        position[i].i = point[i].x - center.x;
        position[i].j = point[i].y - center.y;
        position[i].length = length_v(position[i]);
    }
    /* find shortest and longest distances */
    for(i = 0; i <= size; i++) {
        distance[i] = distance_p(center, point[i]);
        if(distance[i] < shortest) {
            shortest = distance[i];
        }
        if(distance[i] > longest) {
            longest = distance[i];
        }
    }
    /* return the difference in length of the position vectors */
    tao = longest - shortest;
    division_number = (int)sqrt(tao) + 1;
    division_constant = (tao/division_number);
    /* find the points in the current range */
    range[0] = shortest;
    for(i = 0; i < division_number; i++) {
        range[1] = (division_constant * (i + 1)) + shortest;
        /* initialize curve */
        curve[i] = malloc(sizeof(int) * size + 1);
        for(j = 0; j <= size; j++) {
            curve[i][j] = size + 1;
        }
        for(j = 0; j <= size; j++) {
            /* range[0] < distance[j] < range[1] */
            if((distance[j] >= (range[0] - 0.01)) && (distance[j] <= (range[1] + 0.01))) {
                curve[i][j] = j;
            }
        }
        circles[i] = range[0];
        range[0] = range[1];
    }
    circles[i] = range[0];
    fprintf(gnu_files[0], "set parametric\n");
    fprintf(gnu_files[0], "plot %lf*sin(t) + %lf notitle ls %d,", circles[0], center.x, 2);
    fprintf(gnu_files[0], "%lf*cos(t) + %lf notitle ls %d,", circles[0], center.y, 2);
    for(i = 1; i < division_number; i++) {
        fprintf(gnu_files[0], "%lf*sin(t) + %lf notitle ls %d,", circles[i], center.x, 2);
        fprintf(gnu_files[0], "%lf*cos(t) + %lf notitle ls %d,", circles[i], center.y, 2);
    }
    /* printing for debug */
    for(i = 0; i < division_number; i++) {
        printf("curve[%d]: ", i);
        for(j = 0; j <= size; j++) {
            if(curve[i][j] != size + 1) {
                printf("%d ", curve[i][j]);
            }
        }
        printf("\n");
    }//*/
    /* fills neighbors in a 2D array of length (size x 2)
       with the indices of the two neighbors (index of size + 1 means
       a neighbor wasn't found)
       -- runs the contour construction algorithm on a set of points in a curve,
       and stores each point's neighbors
       -- goes through each curve */
    for(i = 0; i < division_number; i++) {
        k = 0;
        /* obtains all points in the current curve */
        for(j = 0; j <= size; j++) {
            if(curve[i][j] != size + 1) {
                contour_point[k] = point[curve[i][j]];
                k++;
            }
        }
        for(j = 0; j < k; j++) {
            printf("contour_point[%d]: %d\n", j, contour_point[j].index);
        }
        //tmp = construct_contour(contour_point, k, gnu_files[2]);
        neighbors = construct_contour(contour_point, k, gnu_files[2]);
        /* finds next k value
        l += k;
        j = 0;
        /* runs from the previous k value to the next k value
        while(m <= l) {
            neighbors[m] = tmp[j];
            printf("neighbors[%d]: %d\n", m, neighbors[m][0], neighbors[m][1]);
            m++;
            j++;
        }
        /* stores previous k value
        m--;*/
    }
    /* plot */
    fprintf(gnu_files[0], "'./gnu_files/lines.tmp' with lines ls 1 title \"neighbors\",");
    fprintf(gnu_files[0], "'./gnu_files/points.tmp' using 1:2 with points pt 7 notitle,");
    fprintf(gnu_files[0], "'' using 1:2:3 with labels point pt 7 offset char -1,-1 notitle\n");
    fprintf(gnu_files[0], "set noparametric\n");
    //free(position);
    //free(distance);
    //free(circles);
    //free(curve);
}

/* constructs triangles given the curve they are on*/
struct triangle_t *construct_triangles(struct point_t **curve, int size,int divisions)
{
    struct triangle_t *triangle = malloc(sizeof(struct triangle_t) * 100);
    return triangle;
}
