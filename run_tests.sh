#!/bin/bash

if (($# < 1)); then
	echo "Usage: ./run_tests.sh [number of random files to test]"
else
	for ((i=1; i <= $1; i++)); do
		./shape_datapoint_generator/random "./datapoints/tests/test$i.dat" "10" &&
		./brute_force "./datapoints/tests/test$i.dat" "./datapoints/tests/test$i.shortest_path"
		./shortest_path "-f" "./datapoints/tests/test$i.dat" "-o" "./datapoints/tests/test$i.calculated_path" &&
		./compare "./datapoints/tests/test$i.calculated_path" "./datapoints/tests/test$i.shortest_path"
	done
fi
