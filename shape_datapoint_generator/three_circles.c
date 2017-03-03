#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct point_t {
    double x;
    double y;
};

struct point_t circle_1(double r, double theta);
struct point_t circle_2(double r, double theta);
struct point_t circle_3(double r, double theta);

int main(int argc, char *argv[])
{
    int size;
    int size_a;
    int size_b;
    FILE *shape = fopen("./datapoints/three_circles.dat", "w");
    if(shape == NULL) {
        printf("File not found. Exiting Program. Good Day.\n");
        exit(EXIT_FAILURE);
    }
    struct point_t *point;
    int i = 0;
    int r = 4;
    if(argc == 1) {
        size = 100;
    }
    else {
        size = atoi(argv[argc - 1]);
    }
    if((size % 3) == 0) {
        size_a = size / 3;
    }
    else {
        size_a = (size / 3) + 1;
    }
    if(size > 100) {
        r *= (size / 10);
    }
    size_b = 2 * size_a;
    point = malloc(sizeof(struct point_t) * size);
    for(; i < size_a; i++) {
        point[i] = circle_1(r, 2*i*M_PI/size_a);
        fprintf(shape, "%d: (%lf, %lf)\n", i, point[i].x, point[i].y);
    }
    for(; i < size_b; i++) {
        point[i] = circle_2(r, 2*(i - size_a)*M_PI/(size_b - size_a));
        fprintf(shape, "%d: (%lf, %lf)\n", i, point[i].x, point[i].y);
    }
    for(; i < size; i++) {
        point[i] = circle_3(r, 2*(i - size_b)*M_PI/(size - size_b));
        fprintf(shape, "%d: (%lf, %lf)\n", i, point[i].x, point[i].y);
    }
    fclose(shape);
    return 0;
}

struct point_t circle_1(double r, double theta)
{
    struct point_t point;
    point.x = r*cos(theta) + r;
    point.y = r*sin(theta) - r;
    return point;
}

struct point_t circle_2(double r, double theta)
{
    struct point_t point;
    point.x = r*cos(theta) - r;
    point.y = r*sin(theta) + r;
    return point;
}

struct point_t circle_3(double r, double theta)
{
    struct point_t point;
    point.x = r*cos(theta) + (2 * r);
    point.y = r*sin(theta) + (2 * r);
    return point;
}
