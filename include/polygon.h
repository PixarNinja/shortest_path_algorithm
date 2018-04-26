/*
 * Definition of the Polygon class and methods used for optimal shape calculations
 *
 * Shortest Path Algorithm
 * Mark Wesley Harris
 * 2018
 */

#ifndef POLYGON_H
#define POLYGON_H

#include <vector>
#include <string>
#include "point.h"

class Polygon {

    private:
        /* calculates the permimeter of the polygon */
        void find_perimeter() {
            int i = 0;
            perimeter = 0.0;

            for(i = 0; i < points.size() - 1; i++) {
                perimeter += distance(points[i], points[i + 1]);
            }
            perimeter += distance(points[i], points[0]);
        }

        /* calculates distance given two points
         * @param P1, the start point
         * @param P2, the end point
         * @return the distance between the points
         */
        double distance(Point P1, Point P2) {
            return sqrt(pow(P2.x - P1.x, 2) + pow(P2.y - P1.y, 2));
        }

        /* calculates the sorted indices of the polygon */
        double find_id() {
            std::vector<Point> tmp = points;
            Point tmp_point;
            int i;
            int j;

            /* bubble sort vertices by index */
            for(i = 0; i < tmp.size(); i++) {
                for(j = tmp.size() - 1; j > i; j--) {
                    if(tmp[j].index > tmp[j - 1].index) {
                        tmp_point = tmp[j];
                        tmp[j] = tmp[j - 1];
                        tmp[j - 1] = tmp_point;
                    }
                }
            }

            id = "";
            for(i = 0; i < tmp.size(); i++) {
                if(i < tmp.size() - 1) {
                    id += std::to_string(tmp[i].index) + "-";
                }
                else {
                    id += std::to_string(tmp[i].index);
                }
            }
        }

    public:
        std::vector<Point> points;
        std::vector<int> shape;
        double perimeter;
        std::string id;

        //////////////////
        // CONSTRUCTORS //
        //////////////////
        
        Polygon();

        /* constructor: base */
        Polygon(std::vector<int> shape, Point *points);

        /* constructor: clone a Polygon */
        Polygon(const Polygon &S);

};

#endif
