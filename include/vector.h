/*
 * Definition of the Vector class and methods used for vector analysis
 *
 * Shortest Path Algorithm
 * Mark Wesley Harris
 * 2018
 */

#pragma once
#include <math.h>
#include "point.h"

class Vector {

    public:
        char *name;
        Point point[2];
        double length;
        double i;
        double j;

        Vector();

        /* constructor */
        Vector(char *name, Point point[2]) {
            this->name = name;
            this->point[0] = point[0];
            this->point[1] = point[1];
            i = 0;
            j = 0;
            find_length();
        }

        /* stores the length of the vector */
        void find_length() {
            length = sqrt(pow(i, 2) + pow(j, 2));
        }
        
        /* finds the angle between another vector
         * @param V, the vector to find the angle with
         * @return the angle in radians
         */
        double find_angle(Vector V) {
            return (acos(dot_product(V) / (length * V.length)));
        }

        /* finds the distance to another vector
         * @param V, the vector to find the distance to
         * @return the distance
         */
        double find_distance(Vector V) {
            return sqrt(pow(V.i - i, 2) + pow(V.j - j, 2));
        }

        /* finds the dot product to another vector
         * @param V, the vector to perform the dot product with
         * @return the dot product
         */
        double dot_product(Vector V) {
            return ((V.i * i) + (V.j * j));
        }

        /* prints information about this vector */
        void print() {
            printf("%s: datapoints = %d (%0.3lf,%0.3lf), %d (%0.3lf,%0.3lf)\n   components = <%0.3lf,%0.3lf>\n   length     = %0.3lf\n\n", name, point[0].index, point[0].x, point[0].y, point[1].index, point[1].x, point[1].y, i, j, length);
        }
        
};
