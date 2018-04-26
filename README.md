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

1. 2000.dat
1. 1000.dat
1. 500.dat
1. 50.dat
1. arrow.dat
1. bigstar.dat
1. bird.dat
1. box.dat
1. butterfly.dat
1. c.dat
1. center.dat
1. circle.dat
1. complex.dat
1. donut.dat
1. ellipse.dat
1. flyingfish.dat
1. frog.dat
1. glob.dat
1. halfbigstar.dat
1. halffrog.dat
1. halfkey.dat
1. halfnested.dat
1. halfseaweed.dat
1. key.dat
1. nested.dat
1. offcenter.dat
1. overlap.dat
1. pill.dat
1. random.dat
1. seaweed.dat
1. simple.dat
1. square.dat
1. star.dat
1. test.dat
1. threecircles.dat
1. tree.dat
1. triangle.dat
1. two_circles.dat

Please note: the GNUplot plotting utility must be installed for the path to be plotted.

## TODO

I need to finish implementing the algorithm; currently, the algorithm calculates
"optimal shapes" within a set of datapoints. I will complete my work on
creating these optimal shapes, and then implement an algorithm that
connects them using Dynamic Programming.
