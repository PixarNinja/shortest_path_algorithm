/*
 * Definition of the Point class and methods used for vector analysis
 *
 * Shortest Path Algorithm
 * Mark Wesley Harris
 * 2018
 */

#ifndef POINT_H
#define POINT_H

#pragma once
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
        Point(double x, double y, int index) {
            this->x = x;
            this->y = y;
            this->index = index;
            theta = 0;
            curvature = 0;
            tao = 0;
            tao_distance = 0;
        }

        /* constructor: clone a Point */
        Point(const Point &P) {
            x = P.x;
            y = P.y;
            theta = P.theta;
            curvature = P.curvature;
            tao = P.tao;
            tao_distance = P.tao_distance;
            index = P.index;
        }

        /* offsets the point
         * @param x, the x offset
         * @param y, the y offset
         */
        void offset(double x_offset, double y_offset) {
            x += x_offset;
            y += y_offset;
        }

        /* tests if a point is equal
         * @param P, the point to test
         * @return true if equal, false otherwise
         */
        bool equals(Point P) {
            if(index == P.index) {
                return true;
            }
            return false;
        }

};

#endif
