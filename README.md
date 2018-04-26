# Shortest Path Algorithm

## Brief

This repo contains the source code for my shortest path algorithm project.
This project involves concepts from vector analysis, computational geometry,
and graph theory which are used to create an algorithm for finding the shortest
path through a set of planar Cartesian coordinates. I plan to fully develop this
algorithm in the hopes of finding a polynomial time solution to the traveling
salesman problem.

This project builds off of the principles of vector analysis I
explored [pixarninja/contour_construction_algorithm](https://github.com/pixarninja/contour_construction_algorithm),
which is an algorithm that generates garunteed shortest paths given that the datapoints
are arranged in a convex shape.

## Purpose

This shortest path algorithm was created to aide in my solution for the 3D
geometric reconstruction of object data from sparse point-clouds. There already exist algorithms
to generate mesh from point-clouds, however the object is not generated correctly
in the case of sparse data. It is my theory that the best possible reconstruction
of any given point-cloud will be the most convex planar mesh that can be generated
from a given contour of the point-cloud; these contours would then be joined together
to create the reconstructed 3D object. My study of this conjecture has lead me to
also theorize that this shape is in fact the shortest halmotonian path for each planar
contour.

## Usage

```
cmake .
make
./shortest_path -f datapoint_file [-o output_file]
```

## Test Files

Use the `run_tests.sh` script to generate random test files
and plot the segments generated from the `shortest_path` and
```brute_force``` executables for comparison. The number of points to generate
for each file is hard-coded at 10 to provide feasible random test cases.

Alternatively, you can enter one of the following files containing
pre-defined datapoint information to run the algorithm on and plot the segments:

0. 2000.dat
0. 1000.dat
0. 500.dat
0. 50.dat
0. arrow.dat
0. bigstar.dat
0. bird.dat
0. box.dat
0. butterfly.dat
0. c.dat
0. center.dat
0. circle.dat
0. complex.dat
0. donut.dat
0. ellipse.dat
0. flyingfish.dat
0. frog.dat
0. glob.dat
0. halfbigstar.dat
0. halffrog.dat
0. halfkey.dat
0. halfnested.dat
0. halfseaweed.dat
0. key.dat
0. nested.dat
0. offcenter.dat
0. overlap.dat
0. pill.dat
0. random.dat
0. seaweed.dat
0. simple.dat
0. square.dat
0. star.dat
0. test.dat
0. threecircles.dat
0. tree.dat
0. triangle.dat
0. two_circles.dat

Please note: the GNUplot plotting utility must be installed for the path to be plotted.

##TODO:

I need to finish implementing the algorithm; currently, the algorithm calculates
"optimal shapes" within a set of datapoints. I will complete my work on
creating these optimal shapes, and then implement an algorithm that
connects them using Dynamic Programming.
