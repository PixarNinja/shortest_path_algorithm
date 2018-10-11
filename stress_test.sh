#!/bin/bash

for ((i=1; i >= 1; i++)); do
	echo "TEST $i:"
	./shape_datapoint_generator/random "./datapoints/tests/stress.dat" "10"
	cat "./datapoints/tests/stress.dat"
	./brute_force "1" "./datapoints/tests/stress.dat"
	./shortest_path "-f" "./datapoints/tests/stress.dat" "-o" "./datapoints/tests/stress.calculated_path"
	./compare "./datapoints/tests/stress.calculated_path" "./datapoints/tests/stress.shortest_path" "./datapoints/tests/stress.dat"
done
