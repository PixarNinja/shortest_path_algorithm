#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct point_t {
    double x;
    double y;
};

struct point_t side_a(double a, double delta);
struct point_t side_b(double b, double delta);
struct point_t side_c(double c, double delta);

int main(int argc, char *argv[])
{
    int size;
    int size_a;
    int size_b;
    FILE *shape = fopen("./datapoints/triangle.dat", "w");
    if(shape == NULL) {
        printf("File not found. Exiting Program. Good Day.\n");
        exit(EXIT_FAILURE);
    }
    struct point_t *point;
    int i = 0;
    double r = -2;
    double s = -2;
    double t = -4;
    if(argc == 1) {
        size = 100;
    }
    else {
        size = atoi(argv[argc - 1]);
    }
    size_a = size / 3;
    size_b = 2 * size_a;
    t += 8.0/size_a;
    point = malloc(sizeof(struct point_t) * size);
    for(; i < size_a; i++) {
        point[i] = side_a(2, r);
        fprintf(shape, "%d: (%lf, %lf)\n", i, point[i].x, point[i].y);
        r += 4.0 / size_a;
    }
    for(; i <= size_b; i++) {
        point[i] = side_b(2, s);
        fprintf(shape, "%d: (%lf, %lf)\n", i, point[i].x, point[i].y);
        s += 4.0 / (size_b - size_a);
    }
    for(; i < size; i++) {
        point[i] = side_c(2, t);
        fprintf(shape, "%d: (%lf, %lf)\n", i, point[i].x, point[i].y);
        t += 8.0 / (size - size_b);
    }
    fclose(shape);
    return 0;
}

struct point_t side_a(double a, double delta)
{
    struct point_t point;
    point.x = delta;
    point.y = point.x - 2;
    return point;
}

struct point_t side_b(double b, double delta)
{
    struct point_t point;
    point.x = delta;
    point.y = 2 - point.x;
    return point;
}

struct point_t side_c(double c, double delta)
{
    struct point_t point;
    point.x = -c;
    point.y = delta;
    return point;
}
