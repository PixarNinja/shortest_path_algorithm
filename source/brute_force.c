#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>

#define SIZE 7

struct point_t {
    double x;
    double y;
    int index;
    double curr;
} point;

void swap(struct point_t *x, struct point_t *y);
void permute(struct point_t *point, struct point_t *shortest, int index, int n);
double distance(struct point_t *point, int first, int last);
int factorial(int n);

int *global_shortest;
int global_count = 0;

int main(int argc, char *argv[])
{
    if(argc != 3) {
        printf("\nUsage: ./shortest [datapoint path] [output path]\n\n");
        exit(EXIT_FAILURE);
    }
    FILE *data = fopen(argv[1], "r");
    FILE *output = fopen(argv[2], "w+");
    FILE *plot;
    struct point_t *point;
    struct point_t *shortest;
    char buf[1024];
    int start = 0;
    int n;
    int size = 0;
    int i = 0;
    const char *gnu_path = "./gnu_files/";
    char *name = malloc(1024);
    char *tmp = malloc(strlen(argv[1]) + 1);
    for(i = 0; i < strlen(argv[1]) + 1; i++)
        tmp[i] = argv[1][i];
    tmp = strtok(tmp, "/");
    while (tmp != NULL) {
        for(i = 0; i < strlen(tmp) + 1; i++)
            name[i] = tmp[i];
        tmp = strtok (NULL, "/");
    }
    name = strtok(name, ".");
    char *path = (char *)malloc(strlen(gnu_path) + strlen(name) + strlen("_shortest_path.gpf") + 1);
    sprintf(path, "%s%s_shortest_path.gpf", gnu_path, name);
    plot = fopen(path, "w+");
    while(fgets(buf, 1024, data)) {
        size++;
    }
    n = size - 1;
    global_shortest = malloc(sizeof(int) * (size + 1));
    for(i = 0; i < size; i++) {
        global_shortest[i] = size;
    }
    i = 0;
    fclose(data);
    data = fopen(argv[1], "r");
    point = malloc(sizeof(struct point_t) * size);
    while(fscanf(data, "%d: (%lf, %lf)", &point[i].index, &point[i].x, &point[i].y) > 0) {
        point[i].curr = DBL_MAX;
        i++;
    }
    shortest = point;
    /* fills array with permutations */
    permute(point, shortest, start, n);
    /* file output */
    printf("\nSHORTEST PATH: %d->", global_shortest[0]);
    fprintf(plot, "%lf %lf\n", point[global_shortest[0]].x, point[global_shortest[0]].y);
    fprintf(output, "%d %d\n", point[global_shortest[0]].index, point[global_shortest[1]].index);
    for(i = 1; i < n; i++) {
        printf("%d->", global_shortest[i]);
        fprintf(plot, "%lf %lf\n", point[global_shortest[i]].x, point[global_shortest[i]].y);
        fprintf(output, "%d %d\n", point[global_shortest[i]].index, point[global_shortest[i + 1]].index);
    }
    printf("%d->%d", global_shortest[i], global_shortest[0]);
    fprintf(plot, "%lf %lf\n", point[global_shortest[i]].x, point[global_shortest[i]].y);
    fprintf(plot, "%lf %lf\n", point[global_shortest[0]].x, point[global_shortest[0]].y);
    fprintf(output, "%d %d\n", point[global_shortest[i]].index, point[global_shortest[0]].index);
    //printf("Total Permutations: %d\n", global_count);
    //printf("Distance: %lf\n", shortest[0].curr);
    //printf("\n");
    fclose(output);
    fclose(data);
    return 0;
}

/* Function to swap values at two pointers */
void swap(struct point_t *x, struct point_t *y)
{
    struct point_t tmp;// = malloc(sizeof(struct point_t));
    tmp.x = x->x;
    tmp.y = x->y;
    tmp.index = x->index;
    x->x = y->x;
    x->y = y->y;
    x->index = y->index;
    y->x = tmp.x;
    y->y = tmp.y;
    y->index = tmp.index;
}

/* calculates all permutations */
void permute(struct point_t *point, struct point_t *shortest, int index, int n)
{
    int i = 0;
    double total = 0;
    double curr = point[0].curr;
    double segment = 0;
    /* base case */
    if(index == n) {
        global_count++;
        /* calculating distance of segments from start to end */
        for(i = 0; i < n; i++) {
            segment = distance(point, i, i + 1);
            total += segment;
        }
        /* calculating the final segment (start and end nodes) */
        segment = distance(point, n, 0);
        total += segment;
        /* checking to see if the most recent total is less than
           the current shortest length */
        if(total < curr) {
            point[0].curr = total;
            global_shortest[0] = shortest[n].index;
            for(i = 0; i <= n; i++) {
                global_shortest[i + 1] = shortest[i].index;
            }
            shortest = point;
            curr = total;
        }
    }
    else {
        for(i = index; i <= n; i++) {
            swap(point + index, point + i);
            permute(point, shortest, index + 1, n);
            swap(point + index, point + i); //backtrack
        }
    }
}

/* calculates distance given index and structure */
double distance(struct point_t *point, int first, int last)
{
    double x_1 = point[first].x;
    double y_1 = point[first].y;
    double x_2 = point[last].x;
    double y_2 = point[last].y;
    double dist = sqrt((x_2 - x_1)*(x_2 - x_1)+(y_2 - y_1)*(y_2 - y_1));
    return dist;
}

/* calculates factorial of an integer */
int factorial(int n)
{
    int result = 1;
    int i = 2;
    for(; i <= n; i++) {
        result *= i;
    }
    return result;
}
