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
    double sign = 0.0;
    /* ensures point is unique */
    while(redo) {
        r_x = (double)(rand() % size);
        r_y = (double)(rand() % size);
        sign = 2 * (double)rand() / (double)RAND_MAX - 1;
        redo = 0;
        for(i = 0; i < size; i++) {
            if(points[i].x == r_x * sign) {
                if(points[i].y == r_y * sign) {
                    redo = 1;
                    break;
                }
            }
        }
    }
    point.x = r_x * sign;
    point.y = r_y * sign;
    return point;
}
