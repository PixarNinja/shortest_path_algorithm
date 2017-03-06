#!/bin/bash

for i in {1..30}
do
	./shape_datapoint_generator/random "./datapoints/tests/test$i.dat" "10" &&
	./tessellate "./datapoints/tests/test$i.dat" "./datapoints/tests/test$i.tes" &&
	./shortest "./datapoints/tests/test$i.dat" "./datapoints/tests/test$i.pth"
	./compare "./datapoints/tests/test$i.tes" "./datapoints/tests/test$i.pth"
done
