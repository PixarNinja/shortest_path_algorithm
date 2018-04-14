/*
 * Definition of the Vector class and methods used for vector analysis
 *
 * Shortest Path Algorithm
 * Mark Wesley Harris
 * 2018
 */

#ifndef VECTOR_H
#define VECTOR_H

#include <math.h>
#include <string>
#include "point.h"

class Vector {

    private:
        /* finds the distance between the endpoints of two vectors
         * @param V1, the first vector
         * @param V2, the second vector
         * @return the distance between the endpoints of the vectors
         */
        double distance_v(Vector V1, Vector V2) {
            return sqrt(pow(V2.i - V1.i, 2) + pow(V2.j - V1.j, 2));
        }

        /* finds the distance between two points
         * @param P1, the first point
         * @param P2, the second point
         * @return the distance between the points
         */
        double distance_p(Point P1, Point P2) {
            return sqrt(pow(P2.x - P1.x, 2) + pow(P2.y - P1.y, 2));
        }

    public:
        std::string name;
        Point start;
        Point end;
        double length;
        double i;
        double j;

        //////////////////
        // CONSTRUCTORS //
        //////////////////
        
        Vector();

        Vector(std::string name, Point start, Point end);

        /* creates a normalized vector centered at start
         * and facing end
         */
        Vector(std::string name, Point start, Point end, int index);

        /////////////
        // METHODS //
        /////////////

        /* refreshes i and j, and stores the length of the vector */
        void refresh();

        /* scoots the vector to start at its endpoint */
        void scoot();

        /* offsets a vector
         * @param x_offset, the x offset
         * @param y_offset, the y offset
         */
        void offset(double x_offset, double y_offset);

        /* prints information about this vector */
        void print();
 
};

#endif
