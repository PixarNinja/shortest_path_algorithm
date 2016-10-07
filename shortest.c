#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#define SIZE 13

struct point_t {
    double x;
    double y;
    char pt;
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
    struct point_t *point = malloc(sizeof(struct point_t) * SIZE);
    struct point_t *shortest = malloc(sizeof(struct point_t) * SIZE);
    shortest = point;
    int start = 0;
    int n = SIZE - 1;
    /* A(0,0) */
    point[0].x = 0;
    point[0].y = 0;
    point[0].pt = 'A';
    point[0].curr = DBL_MAX;
    /* B(1,1) */
    point[1].x = 1;
    point[1].y = 1;
    point[1].pt = 'B';
    point[1].curr = DBL_MAX;
    /* C(-2,1) */
    point[2].x = -2;
    point[2].y = -1;
    point[2].pt = 'C';
    point[2].curr = DBL_MAX;
    /* D(-2,2) */
    point[3].x = -2;
    point[3].y = 2;
    point[3].pt = 'D';
    point[3].curr = DBL_MAX;
    /* E(0,3) */
    point[4].x = 0;
    point[4].y = 3;
    point[4].pt = 'E';
    point[4].curr = DBL_MAX;
    /* F(3,2) */
    point[5].x = 3;
    point[5].y = 2;
    point[5].pt = 'F';
    point[5].curr = DBL_MAX;
    /* G(3,0) */
    point[6].x = 3;
    point[6].y = 0;
    point[6].pt = 'G';
    point[6].curr = DBL_MAX;
    /* H(4,-1) */
    point[7].x = 4;
    point[7].y = -1;
    point[7].pt = 'H';
    point[7].curr = DBL_MAX;
    /* I(1,-3) */
    point[8].x = 1;
    point[8].y = -3;
    point[8].pt = 'I';
    point[8].curr = DBL_MAX;
    /* J(-2,-2) */
    point[9].x = -2;
    point[9].y = -2;
    point[9].pt = 'J';
    point[9].curr = DBL_MAX;
    /* K(-2,-4) */
    point[10].x = -2;
    point[10].y = -4;
    point[10].pt = 'K';
    point[10].curr = DBL_MAX;
    /* L(-4,-1) */
    point[11].x = -4;
    point[11].y = -1;
    point[11].pt = 'L';
    point[11].curr = DBL_MAX;
    /* M(-4,4) */
    point[12].x = -4;
    point[12].y = 4;
    point[12].pt = 'M';
    point[12].curr = DBL_MAX;
    /* fills array with permutations */
    permute(point, shortest, start, n);
    printf("\n");
    printf("Total Permutations: %d\n\n", global_count);
    printf("Shortest Path: %s\n", global_shortest);
    printf("Distance: %lf\n", shortest[0].curr);
    printf("\n");
    return 0;
}

/* Function to swap values at two pointers */
void swap(struct point_t *x, struct point_t *y)
{
    struct point_t tmp;// = malloc(sizeof(struct point_t));
    tmp.x = x->x;
    tmp.y = x->y;
    tmp.pt = x->pt;
    x->x = y->x;
    x->y = y->y;
    x->pt = y->pt;
    y->x = tmp.x;
    y->y = tmp.y;
    y->pt = tmp.pt;
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
            global_shortest[0] = shortest[n].pt;
            for(i = 0; i <= n; i++) {
            	global_shortest[i + 1] = shortest[i].pt;
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
