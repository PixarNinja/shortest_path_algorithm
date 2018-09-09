/*
 * Main program for the shortest path algorithm project
 *
 * Shortest Path Algorithm
 * Mark Wesley Harris
 * 2018
 */

#include <unistd.h>
#include <getopt.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <ctype.h>
#include <unistd.h>

#include <new>
#include <vector>
#include <iostream>
#include <deque>
#include <algorithm>
#include <iterator>
#include <array>
#include <string>

#include "shortest_path.h"

#define NUM_FILES 8

using namespace std;

int main(int argc, char *argv[])
{

    bool print = false;
    FILE *data;
    FILE *output;
    FILE *gnu_files[NUM_FILES];
    vector<Polygon> polygons;
    Polygon tmp_polygon;
    Polygon shortest_path;
    Point *points;
    Point center;
    char buf[1024];
    double range = 0.0;
    double sum_x = 0.0;
    double sum_y = 0.0;
    vector<int *> *segments = new vector<int *> [1];
    vector<int *> edges; //pointer contains the index followed by position
    int **recorded;
    int *mapped;
    int keep_going = 0;
    int size = 0;
    int i = 0;
    int j = 0;
    int k = 0;
    int count = 0;

    /* output file naming */
    char *name = new char[1024];
    char option;
    char *datafile = NULL;
    char *outfile = NULL;

    ////////////
    // GETOPT //
    ////////////

    if(argc == 1) {
        printf("\nUsage: ./shortest_path -f datapoint_path [-p] [-o output_path]\n\n");
        exit(EXIT_FAILURE);
    }

    while ((option = getopt(argc, argv,"hHpf:o:")) != -1) {
        switch (option) {
            case 'f':
            {
                datafile = optarg;
                break;
            }
            case 'o':
            {
                outfile = optarg;
                break;
            }
            case 'p':
            {
                print = true;
                break;
            }
            case 'h':
            case 'H':
            default:
            {
                printf("\nUsage: ./shortest_path -f datapoint_path [-o output_path]\n\n");
                exit(EXIT_FAILURE);
            }
        }
    }

    /* opens file for reading */
    data = fopen(datafile, "r");
    if(!data) {
        printf("\nERROR: unable to open data file. Exiting Program. Good day.\n\n");
        exit(EXIT_FAILURE);
    }
    /* grabs name for output files */
    char *tmp = new char [strlen(datafile) + 1];
    for(i = 0; i < strlen(datafile) + 1; i++)
        tmp[i] = datafile[i];
    tmp = strtok(tmp, "/");
    while (tmp != NULL) {
        for(i = 0; i < strlen(tmp) + 1; i++)
            name[i] = tmp[i];
        tmp = strtok (NULL, "/");
    }
    name = strtok(name, ".");

    /* opens the file for writing */
    if(outfile != NULL) {
        output = fopen(outfile, "w+");
        if(!output) {
            printf("\nERROR: unable to open output file. Exiting Program. Good day.\n\n");
            exit(EXIT_FAILURE);
        }
    }

    /////////////////////
    // OPEN PLOT FILES //
    /////////////////////

    /* setup pathnames for plot output */
    const char *gnu_path = "./gnu_files/";
    char *commands = new char [strlen(gnu_path) + strlen(name) + strlen("_commands.gpf") + 1];
    char *datapoints = new char [strlen(gnu_path) + strlen(name) + strlen("_datapoints.gpf") + 1];
    char *lines = new char [strlen(gnu_path) + strlen(name) + strlen("_lines.gpf") + 1];
    char *extrapoints = new char [strlen(gnu_path) + strlen(name) + strlen("_lines.gpf") + 1];
    char *centerpoint = new char [strlen(gnu_path) + strlen(name) + strlen("_centerpoint.gpf") + 1];
    char *path = new char [strlen(gnu_path) + strlen(name) + strlen("_calculated_path.gpf") + 1];
    char *shortest = new char [strlen(gnu_path) + strlen(name) + strlen("_shortest_path.gpf") + 1];
    char *gnu_tmp = new char [strlen(gnu_path) + strlen("tmp.gpf") + 1];
    sprintf(commands, "%s%s_commands.gpf", gnu_path, name);
    sprintf(datapoints, "%s%s_datapoints.gpf", gnu_path, name);
    sprintf(lines, "%s%s_lines.gpf", gnu_path, name);
    sprintf(extrapoints, "%s%s_extrapoints.gpf", gnu_path, name);
    sprintf(centerpoint, "%s%s_centerpoint.gpf", gnu_path, name);
    sprintf(path, "%s%s_calculated_path.gpf", gnu_path, name);
    sprintf(shortest, "%s%s_shortest_path.gpf", gnu_path, name);
    sprintf(gnu_tmp, "%stmp.gpf", gnu_path);
    gnu_files[0] = fopen(commands, "w+");
    gnu_files[1] = fopen(datapoints, "w+");
    gnu_files[2] = fopen(lines, "w+");
    gnu_files[3] = fopen(extrapoints, "w+");
    gnu_files[4] = fopen(centerpoint, "w+");
    gnu_files[5] = fopen(path, "w+");
    gnu_files[6] = fopen(shortest, "r");
    gnu_files[7] = fopen(gnu_tmp, "w+");
    while(fgets(buf, 1024, data)) {
        size++;
    }
    fclose(data);

    ///////////////////////
    // IMPORT DATAPOINTS //
    ///////////////////////

    points = new Point [size];
    mapped = new int [size];
    recorded = new int * [size];
    /* initializing array */
    for(i = 0; i < size; i++) {
        mapped[i] = 0;
        recorded[i] = new int [size];
        for(j = 0; j < size; j++) {
            recorded[i][j] = 0;
        }
    }
    i = 0;
    data = fopen(datafile, "r");
    while(fscanf(data, "%d: (%lf, %lf)", &points[i].index, &points[i].x, &points[i].y) > 0) {
        if(fabs(points[i].x) > range) {
            range = fabs(points[i].x);
        }
        if(fabs(points[i].y) > range) {
            range = fabs(points[i].y);
        }
        i++;
    }

    /////////////////////
    // PLOT DATAPOINTS //
    /////////////////////

    for(i = 0; i < size; i++) {
        fprintf(gnu_files[1], "%lf %lf %d\n", points[i].x, points[i].y, points[i].index);
    }
    /* plot setup */
    fprintf(gnu_files[0], "set xrange [%lf:%lf]\n", -(range + 1), range + 1);
    fprintf(gnu_files[0], "set yrange [%lf:%lf]\n", -(range + 1), range + 1);
    fprintf(gnu_files[0], "set size ratio 1\n");
    fprintf(gnu_files[0], "set grid\n");
    fprintf(gnu_files[0], "set title \"%s\"\n", datafile);
    fprintf(gnu_files[0], "set style line 1 lc rgb \"black\" lw 1\n");
    fprintf(gnu_files[0], "set style line 2 lc rgb \"red\" lw 3\n");
    fprintf(gnu_files[0], "set style line 3 lc rgb \"#BB0000FF\" lw 6\n");
    fprintf(gnu_files[0], "set style line 4 lc rgb \"#BBFF0000\" lw 6\n");

    // ////////////////////////
    // // CALCULATE POLYGONS //
    // ////////////////////////
    //
    // /* runs experimental algorithm...*/
    // polygons = init_w_polygons(points, size);
    //
    // /* bubble sort polygons by perimeter */
    // for(i = 0; i < polygons.size(); i++) {
    //     for(j = polygons.size() - 1; j > i; j--) {
    //         if(polygons[j].perimeter < polygons[j - 1].perimeter) {
    //             tmp_polygon = polygons[j];
    //             polygons[j] = polygons[j - 1];
    //             polygons[j - 1] = tmp_polygon;
    //         }
    //     }
    // }
    //
    // /* find shortest path */
    // shortest_path = find_shortest_path(polygons, points, size);
    //
    // /* plot shortest path */
    // if(shortest_path.shape.size() > 0) {
    //     printf("\nCALCULATED PATH: ");
    //     for(i = 0; i < shortest_path.shape.size(); i++) {
    //         printf("%d->", points[shortest_path.shape[i]].index);
    //         fprintf(gnu_files[5], "%lf %lf\n", points[shortest_path.shape[i]].x, points[shortest_path.shape[i]].y);
    //     }
    //     printf("%d\n\n", points[shortest_path.shape[0]].index);
    // }
    //
    // ///////////////////////
    // // PLOT POLYGON DATA //
    // ///////////////////////
    //
    // for(Polygon polygon : polygons) {
    //     for(int *segment : polygon.segments) {
    //         if(segment_match(*segments, segment[0], segment[1]) == -1) {
    //             segments->push_back(segment);
    //         }
    //     }
    // }
    //
    // /* print polygon information */
    // printf("\nPOLYGONS:\n");
    // for(i = 0; i < polygons.size(); i++) {
    //     printf("%d: ", i);
    //     for(j = 0; j < (polygons[i]).shape.size(); j++) {
    //         printf("%d ", points[(polygons[i]).shape[j]].index);
    //     }
    //     printf("= %0.2lf\n", polygons[i].perimeter);
    // }

    //////////////////////////
    // CALCULATE W SEGMENTS //
    //////////////////////////

    vector<int *> stored_segments = all_w_segments(points, size);
    for(int *segment : stored_segments) {
        segments->push_back(segment);
    }

    ///////////////////////
    // PLOT SEGMENT DATA //
    ///////////////////////

    /* plot segment information */
    for(i = 0; i < segments->size(); i++) {
        fprintf(gnu_files[2], "%lf %lf %d\n", points[(*segments)[i][0]].x, points[(*segments)[i][0]].y, points[(*segments)[i][0]].index);
        fprintf(gnu_files[2], "%lf %lf %d\n", points[(*segments)[i][1]].x, points[(*segments)[i][1]].y, points[(*segments)[i][1]].index);
        fprintf(gnu_files[2], "\n");
        /* write segments to output file */
        if(outfile != NULL) {
            fprintf(output, "%d %d\n", points[(*segments)[i][0]].index, points[(*segments)[i][1]].index);
        }
    }

    /////////////////////
    // PLOT OTHER DATA //
    /////////////////////

    /* plot perimeter data */
    for(i = 0; i < polygons.size(); i++) {
        sum_x = 0.0;
        sum_y = 0.0;
        for(j = 1; j < (polygons[i]).shape.size(); j++) {
            sum_x += points[(polygons[i]).shape[j]].x;
            sum_y += points[(polygons[i]).shape[j]].y;
        }
        center.x = sum_x / (double)((polygons[i]).shape.size() - 1);
        center.y = sum_y / (double)((polygons[i]).shape.size() - 1);
        //fprintf(gnu_files[3], "%lf %lf %d\n", center.x, center.y, i);
        fprintf(gnu_files[3], "%lf %lf %0.2lf\n", center.x, center.y, polygons[i].perimeter);
    }

    /////////////////////////
    // FINAL PLOT COMMANDS //
    /////////////////////////

    fprintf(gnu_files[0], "plot '%s' using 1:2 with lines ls 1 title \"Tessellations (%ld)\",", lines, segments->size() * 2);
    fprintf(gnu_files[0], "'%s' using 1:2 with points pt 7 notitle,", datapoints);
    fprintf(gnu_files[0], "'' using 1:2:3 with labels point pt 7 offset char -1,-1 notitle,");
    fprintf(gnu_files[0], "'%s' using 1:2:3 with labels point pt 3 offset char -1,-1 notitle, ", extrapoints);
    //fprintf(gnu_files[0], "'%s' using 1:2:3 with labels point pt 2 offset char -1,-1 notitle, ", centerpoint);
    //fprintf(gnu_files[0], "'%s' using 1:2 with lines ls 4 title \"Calculated Path\", ", path);
    fprintf(gnu_files[0], "'%s' using 1:2 with lines ls 3 title \"Shortest Path\"\n", shortest);
    if(outfile != NULL) {
        fclose(output);
    }
    fclose(data);
    fclose(gnu_files[0]);
    fclose(gnu_files[1]);
    fclose(gnu_files[2]);
    fclose(gnu_files[3]);
    fclose(gnu_files[4]);
    fclose(gnu_files[5]);
    fclose(gnu_files[6]);
    char plot[1024];
    sprintf(plot, "gnuplot -persistent %s", commands);
    if(print) {
        system(plot);
    }
    fclose(gnu_files[7]);
    return 0;
}
