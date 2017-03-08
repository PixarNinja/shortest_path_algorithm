#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

struct point_t {
    double x;
    double y;
};

struct point_t generate_random(struct point_t *points, int size);

int main(int argc, char *argv[])
{
    int size;
    if(argc < 3) {
        printf("Usage: [path] [number of datapoints to generate]\n");
        exit(EXIT_FAILURE);
    }
    FILE *shape = fopen(argv[1], "w");
    if(shape == NULL) {
        printf("File not found. Exiting Program. Good Day.\n");
        exit(EXIT_FAILURE);
    }
    struct point_t *points;
    int i = 0;
    /* random seed, only called once */
    srand(time(NULL));
    size = atoi(argv[argc - 1]);
    points = malloc(sizeof(struct point_t) * size);
    for(i = 0; i < size; i++) {
        /* returns a unique randomly generated point */
        points[i] = generate_random(points, size);
        fprintf(shape, "%d: (%lf, %lf)\n", i, points[i].x, points[i].y);
    }
    fclose(shape);
    return 0;
}

struct point_t generate_random(struct point_t *points, int size)
{
    struct point_t point;
    int i = 0;
    int redo = 1;
    double r_x = 0.0;
    double r_y = 0.0;
    double sign_x = 0.0;
    double sign_y = 0.0;
    double epsilon = 0.000001;
    /* ensures point is unique */
    while(redo) {
        r_x = ((double)rand() / (double)RAND_MAX) * size;
        r_y = ((double)rand() / (double)RAND_MAX) * size;
        sign_x = 2 * (double)rand() / (double)RAND_MAX - 1;
        sign_y = 2 * (double)rand() / (double)RAND_MAX - 1;
        if(sign_x >= 0)
            sign_x = 1.0;
        else
            sign_x = -1.0;
        if(sign_y >= 0)
            sign_y = 1.0;
        else
            sign_y = -1.0;
        redo = 0;
        for(i = 0; i < size; i++) {
            if(fabs(points[i].x - (r_x * sign_x)) < epsilon) {
                if(fabs(points[i].y - (r_y * sign_y)) < epsilon) {
                    redo = 1;
                    break;
                }
            }
        }
    }
    point.x = (r_x * sign_x);
    point.y = (r_y * sign_y);
    return point;
}
