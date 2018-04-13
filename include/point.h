/*
 * Definition of the Point class and methods used for vector analysis
 *
 * Shortest Path Algorithm
 * Mark Wesley Harris
 * 2018
 */

#pragma once
#include <math.h>

class Point {

    public:
        double x;
        double y;
        double theta;
        double curvature;
        double tao;
        double tao_distance;
        int index;

        Point();

        /* constructor */
        Point(double x, double y, int index) {
            this->x = x;
            this->y = y;
            this->index = index;
            theta = 0;
            curvature = 0;
            tao = 0;
            tao_distance = 0;
        }

        /* calculates distance given two points
         * @param P, the end point
         */
        double find_distance(struct Point P)
        {
            return sqrt(pow(P.x - x, 2) + pow(P.y - y, 2));
        }

};
