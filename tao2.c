#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#define FILE_NUM 6

struct point_t {
    double x;
    double y;
    int index;
};

struct vector_t {
    struct point_t point[2];
    double length;
    double i;
    double j;
};

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

double shortest_path(struct point_t *point, int size, FILE *gnu_files[FILE_NUM]);
int **construct_neighbors(struct point_t *point, int size, FILE *gnu_files[FILE_NUM]);
struct triangle_t *construct_triangles(struct point_t **curve, int size, int divisions);
double calculate_theta(struct vector_t V1, struct vector_t V2);
double calculate_curvature(struct point_t prev, struct point_t curr, struct point_t next, int size);
double tao_distance(double theta, double length, double curvature);
double distance_p(struct point_t start, struct point_t end);
double distance_v(struct vector_t start, struct vector_t end);
double length_v(struct vector_t v);
struct vector_t subtract_v(struct vector_t V2, struct vector_t V1);
double dot_product(struct vector_t start, struct vector_t end);

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
    distance = shortest_path(point, size - 1, gnu_files);
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
double shortest_path(struct point_t *point, int size, FILE *gnu_files[FILE_NUM])
{
    double total = size;
    int **neighbors = construct_neighbors(point, size, gnu_files);
    return total;
}

/* creates a 2D array (size x 2) of neighbors from the point data */
int **construct_neighbors(struct point_t *point, int size, FILE *gnu_files[FILE_NUM])
{
    struct point_t center;
    struct point_t prev;
    struct vector_t *position = malloc(sizeof(struct vector_t) * size);
    struct vector_t T;
    struct vector_t projection_t;
    struct vector_t projection_n;
    double *distance = malloc(sizeof(double) * size);
    double *circles = malloc(sizeof(double) * size);
    double range[2] = {0};
    double division_constant = 0;
    double shortest = DBL_MAX;
    double longest = 0;
    double tao = 0;
    double sum_x = 0;
    double sum_y = 0;
    int **curve = malloc(sizeof(int *) * size);
    int **neighbors = malloc(sizeof(int *) * size);
    int smallest_curvature[2];
    int i = 0;
    int j = 0;
    int k = 0;
    int index = 0;
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
        curve[i] = malloc(sizeof(int) * size);
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
        k = 0;
    }
    circles[i] = range[0];
    fprintf(gnu_files[0], "set parametric\n");
    fprintf(gnu_files[0], "plot %lf*sin(t) + %lf notitle ls %d,", circles[0], center.x, 2);
    fprintf(gnu_files[0], "%lf*cos(t) + %lf notitle ls %d,", circles[0], center.y, 2);
    fprintf(gnu_files[0], "%lf*sin(t) + %lf notitle ls %d,", circles[1], center.x, 2);
    fprintf(gnu_files[0], "%lf*cos(t) + %lf notitle ls %d,", circles[1], center.y, 2);
    fprintf(gnu_files[0], "%lf*sin(t) + %lf notitle ls %d,", circles[2], center.x, 2);
    fprintf(gnu_files[0], "%lf*cos(t) + %lf notitle ls %d,", circles[2], center.y, 2);
    fprintf(gnu_files[0], "%lf*sin(t) + %lf notitle ls %d,", circles[3], center.x, 2);
    fprintf(gnu_files[0], "%lf*cos(t) + %lf notitle ls %d,", circles[3], center.y, 2);
    /* printing for debug
    for(i = 0; i < division_number; i++) {
        printf("curve[%d]: ", i);
        for(j = 0; j <= size; j++) {
            if(curve[i][j] != size + 1) {
                printf("%d ", curve[i][j]);
            }
        }
        printf("\n");
    }*/
    /* fills neighbors in a 2D array of length (size x 2)
       with the indices of the two neighbors (index of size + 1 means
       a neighbor wasn't found)
       -- initializes neighbors[] array to a non-existant node */
    for(i = 0; i <= size; i++) {
        /* each point can have 0, 1, or 2 neighbors */
        neighbors[i] = malloc(sizeof(int) * 2);
        neighbors[i][0] = size + 1;
        neighbors[i][1] = size + 1;
    }
    for(i = 0; i < division_number; i++) {
        /* calculate distances between all paths within curve[i] */
        for(j = 0; j <= size; j++) {
            /* initializes smallest_curvature[] */
            smallest_curvature[0] = size + 1;
            smallest_curvature[1] = size + 1;
            /* initializes T */
            T.i = position[j].j;
            T.j = -position[j].i;
            T.length = length_v(T);
            T.point[0].x = center.x;
            T.point[0].y = center.y;
            T.point[0].index = size + 1;
            T.point[1].x = center.x + T.i;
            T.point[1].y = center.y + T.j;
            T.point[1].index = size + 1;
            /* store vectors_t data */
            fprintf(gnu_files[3], "%lf %lf %lf %lf\n", center.x, center.y, T.i, T.j);
            /* store vectors_n data */
            fprintf(gnu_files[4], "%lf %lf %lf %lf\n", center.x, center.y, position[j].i, position[j].j);
            /* first find a starting point (index j) */
            if(curve[i][j] != (size + 1)) {
                /* initializing previous point for j */
                prev.x = point[j].x + (position[j].i / position[j].length);
                prev.y = point[j].y + (position[j].j / position[j].length);
                prev.index = size + 1;
                /* then compare distances between all the other
                   points (index k) */
                for(k = 0; k <= size; k++) {
                    if((curve[i][k] != (size + 1)) && (k != j)) {
                        /* stores the two shortest paths
                           -- stores paths if the current curvature is
                              less than the curvature of the last
                              point */
                        /* vector from point k to T endpoint */
                        projection_t.point[0].x = position[k].point[1].x;
                        projection_t.point[0].y = position[k].point[1].y;
                        projection_t.point[0].index = position[k].point[1].index;
                        projection_t.point[1].x = T.point[1].x;
                        projection_t.point[1].y = T.point[1].y;
                        projection_t.point[1].index = T.point[1].index;
                        projection_t.i = projection_t.point[1].x - projection_t.point[0].x;
                        projection_t.j = projection_t.point[1].y - projection_t.point[0].y;
                        projection_t.length = length_v(projection_t);
                        /* vector from point k to N endpoint */
                        projection_n.point[0].x = position[k].point[1].x;
                        projection_n.point[0].y = position[k].point[1].y;
                        projection_n.point[1].x = position[j].point[1].x;
                        projection_n.point[1].y = position[j].point[1].y;
                        projection_n.i = projection_n.point[1].x - projection_n.point[0].x;
                        projection_n.j = projection_n.point[1].y - projection_n.point[0].y;
                        projection_n.length = length_v(projection_n);
                        /* check if projection_t is the edge case */
                        printf("i = %d, j = %d, k = %d\n", i, j, k);
                        printf("projection_t.length: %lf, sqrt(2 * T.length) = %lf\n", projection_t.length, sqrt(2 * T.length));
                        if(isgreaterequal(projection_t.length, sqrt(2 * T.length)) && islessequal(projection_t.length, sqrt(2 * T.length))) {
                            /* N-axis
                               -- check if projection_n is the edge case */
                            if(isgreaterequal(projection_n.length, sqrt(2 * position[j].length)) && isgreaterequal(projection_n.length, sqrt(2 * position[j].length))) {
                                printf("ERROR: duplicated point detected. Exiting Program.\n");
                                exit(EXIT_FAILURE);
                            }
                            /* store in "negative" array */
                            else if(isgreaterequal(projection_n.length, sqrt(2 * position[j].length))) {
                                index = 1;
                            }
                            /* store in "positive" array */
                            else {
                                index = 0;
                            }
                        }
                        /* T-axis
                           -- store in "negative" array */
                        else if(isgreaterequal(projection_t.length, sqrt(2 * T.length))) {
                            index = 1;
                        }
                        else {
                            /* store in "positive" array */
                            index = 0;
                        }
                        /* curvature calculations */
                        if(smallest_curvature[index] == (size + 1)) {
                            smallest_curvature[index] = k;
                        }
                        else if (islessequal(calculate_curvature(prev, point[j], point[k], size), calculate_curvature(prev, point[j], point[smallest_curvature[index]], size))) {
                            smallest_curvature[index] = k;
                        }
                    }
                }
                /* printing for debug
                printf("curv[%d] = %0.2lf, : curv[%d] = %0.2lf\n", smallest_curvature[0], calculate_curvature(prev, point[j], point[smallest_curvature[0]]), smallest_curvature[1], calculate_curvature(prev, point[j], point[smallest_curvature[1]])); */
                /* sets neighbors of node i */
                neighbors[j][0] = smallest_curvature[0];
                neighbors[j][1] = smallest_curvature[1];
            }
        }
    }
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
    /* plot */
    fprintf(gnu_files[0], "'./gnu_files/vectors_n.tmp' with vectors arrowstyle 1 title \"N vectors\",");
    fprintf(gnu_files[0], "'./gnu_files/vectors_t.tmp' with vectors arrowstyle 2 title \"T vectors\",");
    fprintf(gnu_files[0], "'./gnu_files/lines.tmp' with lines ls 1 title \"neighbors\",");
    fprintf(gnu_files[0], "'./gnu_files/points.tmp' using 1:2 with points pt 7 notitle,");
    fprintf(gnu_files[0], "'' using 1:2:3 with labels point pt 7 offset char -1,-1 notitle\n");
    fprintf(gnu_files[0], "set noparametric\n");
    //free(position);
    //free(distance);
    //free(circles);
    //free(curve);
    //free(neighbors);
    return curve;
}

/* constructs triangles given the curve they are on*/
struct triangle_t *construct_triangles(struct point_t **curve, int size,int divisions)
{
    struct triangle_t *triangle = malloc(sizeof(struct triangle_t) * 100);
    return triangle;
}

/* calculates curvature given the previous, current, and next points */
double calculate_curvature(struct point_t prev, struct point_t curr, struct point_t next, int size)
{
    /* initializing point data
       -- initializing vector T1 */
    struct vector_t T1;
    T1.point[0].x = curr.x;
    T1.point[0].y = curr.y;
    T1.point[0].index = curr.index;
    T1.point[1].x = prev.x;
    T1.point[1].y = prev.y;
    T1.point[1].index = prev.index;
    T1.i = T1.point[1].x - T1.point[0].x;
    T1.j = T1.point[1].y - T1.point[0].y;
    T1.length = length_v(T1);
    /* initializing vector V2, <curr, next> */
    struct vector_t V2;
    V2.point[0].x = curr.x;
    V2.point[0].y = curr.y;
    V2.point[0].index = curr.index;
    V2.point[1].x = next.x;
    V2.point[1].y = next.y;
    V2.point[1].index = next.index;
    V2.i = V2.point[1].x - V2.point[0].x;
    V2.j = V2.point[1].y - V2.point[0].y;
    V2.length = length_v(V2);
    /* -- initializing vector T2 */
    struct vector_t T2;
    T2.point[0].x = curr.x;
    T2.point[0].y = curr.y;
    T2.point[0].index = curr.index;
    T2.point[1].x = (V2.i / V2.length) + T2.point[0].x;
    T2.point[1].y = (V2.j / V2.length) + T2.point[0].y;
    T2.point[1].index = size + 1;
    T2.i = T2.point[1].x - T2.point[0].x;
    T2.j = T2.point[1].y - T2.point[0].y;
    T2.length = length_v(T2);
    return (calculate_theta(T1, T2));
}

/* calculates theta given two vectors */
double calculate_theta(struct vector_t V1, struct vector_t V2)
{
    return acos(dot_product(V1, V2)/(V1.length * V2.length));
}

/* calculates dot product of two vectors */
double dot_product(struct vector_t start, struct vector_t end)
{
    double i_1 = start.i;
    double j_1 = start.j;
    double i_2 = end.i;
    double j_2 = end.j;
    return ((i_2 * i_1) + (j_2 * j_1));
}

/* calculates tao_distance given theta, length and curvature */
double tao_distance(double theta, double length, double curvature)
{
    return abs((tan(theta) + 1) * length / curvature);
}

/* calculates distance given two points */
double distance_p(struct point_t start, struct point_t end)
{
    double x_1 = start.x;
    double y_1 = start.y;
    double x_2 = end.x;
    double y_2 = end.y;
    return sqrt((x_2 - x_1) * (x_2 - x_1) + (y_2 - y_1) * (y_2 - y_1));
}

/* calculates distance given two vectors */
double distance_v(struct vector_t start, struct vector_t end)
{
    double i_1 = start.i;
    double j_1 = start.j;
    double i_2 = end.i;
    double j_2 = end.j;
    return sqrt((i_2 - i_1) * (i_2 - i_1) + (j_2 - j_1) * (j_2 - j_1));
}

/* calculates length of a single vector */
double length_v(struct vector_t V)
{
    double i = V.i;
    double j = V.j;
    return sqrt((i * i) + (j * j));
}

/* calculates difference of a two vectors */
struct vector_t subtract_v(struct vector_t V2, struct vector_t V1)
{
    struct vector_t S;
    S.point[0].x = V1.point[1].x;
    S.point[0].y = V1.point[1].y;
    S.point[1].x = V2.point[1].x;
    S.point[1].y = V2.point[1].y;
    S.i = S.point[1].x - S.point[0].x;
    S.j = S.point[1].y - S.point[0].y;
    S.length = length_v(S);
    return S;
}
