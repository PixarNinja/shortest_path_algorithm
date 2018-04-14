/*
 * Definition of the Point class and methods used for vector analysis
 *
 * Shortest Path Algorithm
 * Mark Wesley Harris
 * 2018
 */

#ifndef POINT_H
#define POINT_H

#include <math.h>

class Point {

    private:
        /* calculates distance given two points
         * @param P1, the start point
         * @param P2, the end point
         * @return the distance between the points
         */
        double distance(Point P1, Point P2) {
            return sqrt(pow(P2.x - P1.x, 2) + pow(P2.y - P1.y, 2));
        }

    public:
        double x;
        double y;
        double theta;
        double curvature;
        double tao;
        double tao_distance;
        int index;

        Point();

        /* constructor: base */
        Point(double x, double y, int index);

        /* constructor: verbose */
        Point(double x, double y, int index, double theta, double curvature, double tao, double tao_distance);

        /* constructor: clone a Point */
        Point(const Point &P);

        /* offsets the point
         * @param x, the x offset
         * @param y, the y offset
         */
        void offset(double x_offset, double y_offset);

        /* tests if a point is equal
         * @param P, the point to test
         * @return true if equal, false otherwise
         */
        bool equals(Point P);

};

#endif
