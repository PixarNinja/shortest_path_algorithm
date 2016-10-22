#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#define SIZE 13

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

char global_shortest[SIZE + 2] = {'\0'};
int global_count = 0;

int main(void)
{
    FILE *data = fopen("./datapoints/datapoints.dat", "r");
    struct point_t *point;
    struct point_t *shortest = malloc(sizeof(struct point_t) * SIZE);
    char buf[1024];
    int start = 0;
    int n = SIZE - 1;
    int size = 0;
    int i = 0;
    while(fgets(buf, 1024, data)) {
        size++;
    }
    fclose(data);
    data = fopen("./datapoints/datapoints.dat", "r");
    point = malloc(sizeof(struct point_t) * size);
    while(fscanf(data, "%d: (%lf, %lf)", &point[i].index, &point[i].x, &point[i].y) > 0) {
        point[i].curr = DBL_MAX;
        i++;
    }
    fclose(data);
    shortest = point;
    /* fills array with permutations */
    permute(point, shortest, start, n);
    printf("\n");
    printf("Total Permutations: %d\n\n", global_count);
    printf("Shortest Path: %s\n", global_shortest);
    printf("Distance: %lf\n", shortest[0].curr);
    printf("\n");
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
            /* storing new shortest path */
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
            swap(point+index, point+i);
            permute(point, shortest, index+1, n);
            swap(point+index, point+i); //backtrack
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
