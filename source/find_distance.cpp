#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include "point.h"

using namespace std;

/* calculates distance given two points
 * @param P1, the start point
 * @param P2, the end point
 * @return the distance between the points
 */
double distance(Point P1, Point P2) {
    return sqrt(pow(P2.x - P1.x, 2) + pow(P2.y - P1.y, 2));
}

int main(int argc, char *argv[])
{
    if(argc != 3) {
        printf("\nUsage: ./find_distance [datapoint path] [path]\n");
        printf("E.g. './find_distance my_path.dat  1,2,3,4,5'\n\n");
        exit(EXIT_FAILURE);
    }
    FILE *data = fopen(argv[1], "r");
    Point *points;
    char buf[1024];
    int start = 0;
    int size = 0;
    int i = 0;
    while(fgets(buf, 1024, data)) {
        size++;
    }
    fclose(data);
    data = fopen(argv[1], "r");
    points = new Point [size];
    while(fscanf(data, "%d: (%lf, %lf)", &points[i].index, &points[i].x, &points[i].y) > 0) {
        i++;
    }

    string s = string(argv[2]);
    string delimiter = ",";
    int *order = new int [size];

    size_t pos = 0;
    std::string token;
    i = 0;
    while ((pos = s.find(delimiter)) != std::string::npos) {
        token = s.substr(0, pos);
        order[i++] = stoi(token);
        s.erase(0, pos + delimiter.length());
    }

    double sum = 0.0;
    for(i = 1; i < size; i++) {
        sum += distance(points[order[i - 1]], points[order[i]]);
    }

    printf("\nDISTANCE: %lf\n\n", sum);

    return 0;
}
