#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct point_t {
    double x;
    double y;
};

struct point_t ellipse(double a, double b, double theta);

int main(int argc, char *argv[])
{
    int size;
    FILE *shape = fopen("./datapoints/ellipse.dat", "w");
    if(shape == NULL) {
        printf("File not found. Exiting Program. Good Day.\n");
        exit(EXIT_FAILURE);
    }
    struct point_t *point;
    int i = 0;
    if(argc == 1) {
        size = 100;
    }
    else {
        size = atoi(argv[argc - 1]);
    }
    point = malloc(sizeof(struct point_t) * size);
    for(; i < size; i++) {
        point[i] = ellipse(4, 2, 2*i*M_PI/size);
        fprintf(shape, "%d: (%lf, %lf)\n", i, point[i].x, point[i].y);
    }
    fclose(shape);
    return 0;
}

struct point_t ellipse(double a, double b, double theta)
{
    struct point_t point;
    point.x = a*cos(theta);
    point.y = b*sin(theta);
    return point;
}
