#include "construct_contour.h"

/* calculates the contour path and returns it in an array
 * of neighbors for each point */
int **construct_contour(struct point_t *points, int size, FILE *plot)
{
    struct vector_t initial;
    struct vector_t check;
    struct point_t best;
    struct point_t start = points[0];
    struct point_t end = start;
    struct point_t prev = start;
    struct point_t center;
    double total = 0.0;
    double sum_x = 0.0;
    double sum_y = 0.0;
    double theta = DBL_MAX;
    double tmp = 0.0;
    int **neighbors = calloc(size + 1, sizeof(int *) * (size + 1));
    int global_count = 0;
    int n = 0;
    int i = 0;
    int j = 0;
    int index = 0;
    int count = 0;
    int total_size = size - 1;
    int num_points;
    for(i = 0; i <= size; i++) {
        neighbors[i] = malloc(sizeof(int) * 2);
        neighbors[i][0] = size + 1;
        neighbors[i][1] = size + 1;
    }
    best.x = DBL_MAX;
    best.y = DBL_MAX;
    best.tao_d = DBL_MAX;
    best.index = INT_MAX;
    /* calculate average point */
    for(i = 0; i <= size; i++) {
        sum_x += points[i].x;
        sum_y += points[i].y;
    }
    center.x = sum_x / size;
    center.y = sum_y / size;
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
    /* initializing structure curr */
    struct point_t *curr = malloc(sizeof(struct point_t) * (size + 1));
    for(i = 0; i <= size; i++) {
        curr[i].x = DBL_MAX;
        curr[i].y = DBL_MAX;
        curr[i].tao_d = DBL_MAX;
        curr[i].index = INT_MAX;
    }
    struct curvature_t k;
    /* initializing point 1 data for structure k
    -- initializing vector V */
    k.V.point[0].x = start.x;
    k.V.point[0].y = start.y;
    k.V.point[0].index = start.index;
    k.V.i = 0;
    k.V.j = 0;
    k.V.length = 0;
    /* -- initializing vector T1 */
    k.T1.point[0].x = start.x;
    k.T1.point[0].y = start.y;
    k.T1.point[0].index = start.index;
    /* checks for the point with the smallest deviation in angle from
       the position vector of the starting point (selects the "second"
       point) */
    for(i = 0; i <= size; i ++) {
        /* skip the starting point */
        if(points[i].index == start.index) {
            continue;
        }
        check.point[0].x = start.x;
        check.point[0].y = start.y;
        check.point[0].index = start.index;
        check.point[1].x = points[i].x;
        check.point[1].y = points[i].y;
        check.point[1].index = points[i].index;
        check.i = check.point[1].x - check.point[0].x;
        check.j = check.point[1].y - check.point[0].y;
        check.length = length_v(check);
        tmp = angle_v(initial, check);
        if(tmp < theta) {
            theta = tmp;
            index = i;
        }
    }
    /* points T1 in the direction of the "second point" */
    k.T1.i = (points[index].x - start.x) / distance_p(points[index], start);
    k.T1.j = (points[index].y - start.y) / distance_p(points[index], start);
    k.T1.point[1].x = start.x + k.T1.i;
    k.T1.point[1].y = start.y + k.T1.j;
    k.T1.point[1].index = INT_MAX;
    k.T1.length = 1;
    /* -- initializing vector T2 */
    k.T2.point[0].x = start.x;
    k.T2.point[0].y = start.y;
    k.T2.point[0].index = start.index;
    index = 0;
    /* plot */
    fprintf(plot, "%lf %lf %d\n", start.x, start.y, start.index);
    /* remove the start point from the list */
    for(i = 0; i <= size; i++) {
        points[i] = points[i + 1];
    }
    size--;
    /* outer loop, calculates total distance */
    while(global_count <= total_size) {
        i = 0;
        /* refreshing best index */
        best.tao_d = DBL_MAX;
        best.index = start.index;
        /* loops through all possible indices from start */
        while(count < size) {
            /* initializing point 2 data for structure k
               -- initializing vector V */
            k.V.point[1].x = points[i].x;
            k.V.point[1].y = points[i].y;
            k.V.point[1].index = points[i].index;
            k.V.i = k.V.point[1].x - k.V.point[0].x;
            k.V.j = k.V.point[1].y - k.V.point[0].y;
            k.V.length = length_v(k.V);
            /* -- initializing vector T2 */
            k.T2.point[1].x = k.V.point[1].x;
            k.T2.point[1].y = k.V.point[1].y;
            k.T2.point[1].index = INT_MAX;
            k.T2.i = (k.T2.point[1].x - k.T2.point[0].x) / k.V.length;
            k.T2.j = (k.T2.point[1].y - k.T2.point[0].y) / k.V.length;
            k.T2.point[1].x = k.V.point[0].x + k.T2.i;
            k.T2.point[1].y = k.V.point[0].y + k.T2.j;
            k.T2.length = length_v(k.T2);
            /* -- initializing tao, theta, and curvature */
            k.tao = (dot_product(k.T1, k.T2)); //length of T1 and T2 is always 1
            if(k.tao <= -1.0) {
                k.tao = -1.0;
            }
            else if(k.tao >= 1.0) {
                k.tao = 1.0;
            }
            k.theta = calculate_theta(k);
            k.curvature = calculate_curvature(k);
            /* initializing structure curr */
            curr[i].x = k.V.point[1].x;
            curr[i].y = k.V.point[1].y;
            curr[i].index = k.V.point[1].index;
            /* calculating tao-distance */
            curr[i].tao_d = tao_distance(k);
            k.V.point[1].tao_d = curr[i].tao_d;
            i++;
            count++;
        }
        neighbors[index][0] = prev.index;
        /* sets the previous point as the previous best point */
        prev = best;
        /* find point with the lowest tao-distance */
        for(i = 0; i <= size; i++) {
            if(best.tao_d > curr[i].tao_d) {
                best.x = curr[i].x;
                best.y = curr[i].y;
                best.index = curr[i].index;
                best.tao_d = curr[i].tao_d;
                n = i;
            }
        }
        /* remove the best point from the list */
        for(i = 0; i <= size; i++) {
            if(curr[i].index == best.index) {
                for(j = i; j <= size; j++) {
                    points[j] = points[j + 1];
                }
                size--;
                break;
            }
        }
        neighbors[index][1] = best.index;
        total += distance_p(start, best);
        /* plot */
        fprintf(plot, "%lf %lf %d\n", best.x, best.y, best.index);
        /* reinitializing structure k
           -- reinitializing vector V */
        k.V.point[1].x = best.x;
        k.V.point[1].y = best.y;
        k.V.point[1].index = best.index;
        k.V.i = k.V.point[1].x - k.V.point[0].x;
        k.V.j = k.V.point[1].y - k.V.point[0].y;
        k.V.length = length_v(k.V);
        /* -- reinitializing vector T1 */
        k.T2.point[1].x = best.x;
        k.T2.point[1].y = best.y;
        k.T2.point[1].index = INT_MAX;
        k.T2.i = (k.T2.point[1].x - k.T2.point[0].x) / k.V.length;
        k.T2.j = (k.T2.point[1].y - k.T2.point[0].y) / k.V.length;
        k.T2.length = length_v(k.T2);
        k.T1.point[0].x = best.x;
        k.T1.point[0].y = best.y;
        k.T1.point[0].index = best.index;
        k.T1.point[1].x = best.x + k.T2.i;
        k.T1.point[1].y = best.y + k.T2.j;
        k.T1.point[1].index = INT_MAX;
        k.T1.i = (k.T1.point[1].x - k.T1.point[0].x);
        k.T1.j = (k.T1.point[1].y - k.T1.point[0].y);
        k.T1.length = length_v(k.T1);
        /* -- shift starting point to best point */
        start = best;
        /* -- initializing vector T2 */
        k.T2.point[0].x = start.x;
        k.T2.point[0].y = start.y;
        k.T2.point[0].index = start.index;
        k.T2.i = 0;
        k.T2.j = 0;
        k.T2.length = 0;
        /* -- initializing vector V */
        k.V.point[0].x = start.x;
        k.V.point[0].y = start.y;
        k.V.point[0].index = start.index;
        k.V.i = 0;
        k.V.j = 0;
        k.V.length = 0;
        count = 0;
        global_count++;
        index++;
    }
    neighbors[0][0] = best.index;
    neighbors[0][0] = best.index;
    /* final point */
    fprintf(plot, "%lf %lf %d\n", end.x, end.y, end.index);
    fprintf(plot, "\n");
    total += distance_p(best, end);
    return neighbors;
}

/* calculates curvature given structure k */
double calculate_curvature(struct curvature_t k)
{
    return (sqrt(pow((k.T2.i - k.T1.i), 2) + pow((k.T2.j - k.T1.j), 2)) / calculate_theta(k));
}

/* calculates theta given structure k */
double calculate_theta(struct curvature_t k)
{
    return (acos(k.tao) + (M_PI / 180));
}

/* calculates distance given index and structure */
double tao_distance(struct curvature_t k)
{
    return (k.V.length + k.curvature);
}

/* calculates angle between two vectors */
double angle_v(struct vector_t start, struct vector_t end)
{
    return (acos(dot_product(start, end) / (start.length * end.length)));
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

/* calculates length of a single vectors */
double length_v(struct vector_t v)
{
    double i = v.i;
    double j = v.j;
    return sqrt((i * i) + (j * j));
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

/* prints structure k for debugging */
void print_k(struct curvature_t k)
{
    printf("V: %d[0](%lf, %lf), %d[1](%lf, %lf), <%lf, %lf>, |V| = %lf\n", k.V.point[0].index, k.V.point[0].x, k.V.point[0].y, k.V.point[1].index, k.V.point[1].x, k.V.point[1].y, k.V.i, k.V.j, k.V.length);
    printf("T1: point[0](%lf, %lf), point[1](%lf, %lf), <%lf, %lf>, |T1| = %lf\n", k.T1.point[0].x, k.T1.point[0].y, k.T1.point[1].x, k.T1.point[1].y, k.T1.i, k.T1.j, k.T1.length);
    printf("T2: point[0](%lf, %lf), point[1](%lf, %lf), <%lf, %lf>, |T2| = %lf\n", k.T2.point[0].x, k.T2.point[0].y, k.T2.point[1].x, k.T2.point[1].y, k.T2.i, k.T2.j, k.T2.length);
    printf("curvature: %lf; ", k.curvature);
    printf("angle = %lf; ", k.theta * 180 / M_PI);
    printf("tao = %lf; ", k.tao);
    printf("tao-distance = %lf\n\n", k.V.point[1].tao_d);
}
