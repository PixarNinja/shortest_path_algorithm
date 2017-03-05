#!/bin/bash

./shape_datapoint_generator/random ./datapoints/random.dat 10 && ./tessellate datapoints/random.dat && ./shortest ./datapoints/random.dat
