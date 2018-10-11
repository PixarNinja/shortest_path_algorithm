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
#include "vector.h"

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
                    else if(tmp[j].index == tmp[j - 1].index) {
                        tmp.erase(tmp.begin() + j);
                        i = 0;
                        j = 1;
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

        /* returns the determinant between two vectors
         * @param V1, the first vector
         * @param V2, the second vector
         * @return the determinant, i_1*j_2 - j_1*i_2
         */
        double determinant(Vector V1, Vector V2) {
            return ((V1.i * V2.j) - (V2.i * V1.j));
        }

        /* finds the angle between another vector
         * @param V1, the first vector
         * @param V2, the second vector
         * @return the angle in radians
         */
        double angle(Vector V1, Vector V2) {
            return (acos(dot_product(V1, V2) / (V1.length * V2.length)));
        }

        /* finds the dot product of two vectors
         * @param V1, the first vector
         * @param V2, the second vector
         * @return the dot product of the vectors
         */
        double dot_product(Vector V1, Vector V2) {
            return ((V1.i * V2.i) + (V1.j * V2.j));
        }

        /* finds the projection_V1(V2) of two vectors
         * @param V1, the first vector
         * @param V2, the second vector
         * @return projection_V1(V2)
         */
        Vector projection(Vector V1, Vector V2) {
            double coeff = (dot_product(V1, V2) / (V1.length * V1.length));
            double i = V1.i * coeff;
            double j = V1.j * coeff;
            Vector result = Vector(V1);
            result.end = Point(result.start);
            result.end.offset(i, j);
            result.refresh();
            return result;
        }

    public:
        std::vector<Point> points;
        std::vector<int> shape;
        std::vector<int *> segments;
        double perimeter;
        std::string id;

        //////////////////
        // CONSTRUCTORS //
        //////////////////
        
        Polygon();

        /* constructor: base */
        Polygon(std::vector<int> shape, Point *points);

        /* constructor: initialize with segments */
        Polygon(std::vector<int> shape, std::vector<int *> segments, Point *points);

        /* constructor: clone a Polygon */
        Polygon(const Polygon &S);

        /////////////
        // METHODS //
        /////////////

        /* tests if a point is part of a point array
         * @param points, the point array to test
         * @param size, the size of the array
         * @param vertex, the vertex value to find
         * @return the index of the point if found, -1 otherwise
         */
        int point_match(Point *points, int size, int vertex);

        /* tests if a segment is part of a segement vector
         * @param segements, the segment vector to test
         * @param beginning, the beginning index
         * @param end, the end index
         * @return the index of the segment if found, -1 otherwise
         */
        int segment_match(std::vector<int *> segments, int beginning, int end);

        /* re-orders the points of the shape into a hull
         * @param data, the overall point array
         * @param size, the size of the array
         */
        void create_hull(Point *data, int size);

        /* finds if a point is inside of the shape (requires the shape to be sorted as a hull)
         * @param test, the point to test
         * @param points, the point array in which the test point exists
         * @return true if the point is inside the shape, false otherwise
         */
        bool contains(Point test, Point *points, int size);

        /* tests if a polygon is equal
         * @param S, the polygon to test
         * @return true if equal, false otherwise
         */
        bool equals(Polygon S);

};

#endif
