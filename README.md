#Shortest Path Algorithm

##Usage:

```
make
./tesselate [datapoint file]
```
Enter one of the following files containing datapoint information
to run the algorithm on it:

0. arrow.dat
0. bigstar.dat
0. circle.dat
0. complex.dat
0. donut.dat
0. ellipse.dat
0. frog.dat
0. nested.dat
0. offcenter.dat
0. pill.dat
0. simple.dat
0. square.dat
0. star.dat
0. three_circles.dat
0. tree.dat
0. triangle.dat
0. two_circles.dat

Please note: the GNUplot plotting utility must be installed for the path to be plotted.

##Brief:

This algorithm calculates the shortest path between a set of data points
using the Contour Construction Algorithm (https://github.com/pixarninja/contour_construction_algorithm) and principles
of geometry.

##TODO:

Finish implementing the algorithm; currently, the algorithm calculates
"optimal shapes" within a set of datapoints. I will complete my work on
creating these optimal shapes, and then implement an algorithm that
connects them using Dynamic Programming.

Note that the algorithm works for all "round" datasets, since it follows
directly from my work with the Contour Construction Algorithm. Thus
the Best Case Time Complexity of this algorithm is O(n<sup>2</sup>).
