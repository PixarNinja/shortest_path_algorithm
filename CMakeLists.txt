cmake_minimum_required(VERSION 3.0.2)
project (shortest_path_algorithm)

set( CMAKE_EXPORT_COMPILE_COMMANDS on )

include_directories(include)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --std=c++11 -g -Wpedantic")
include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_executable(brute_force source/brute_force.c)
add_executable(compare source/compare.c)
#add_executable(multiple source/multiple_contours.cpp)
add_executable(plot source/plot.c)
#add_executable(pop source/pop_reverse.cpp)
#add_executable(tessellate source/tessellate.cpp)
add_executable(find_distance source/find_distance.cpp source/point.cpp)
add_executable(shortest_path source/main.cpp source/shortest_path.cpp source/vector.cpp source/point.cpp)

target_link_libraries(brute_force m)
#target_link_libraries(multiple m)
#target_link_libraries(pop m)
#target_link_libraries(tessellate m)
target_link_libraries(find_distance m)
target_link_libraries(shortest_path m)
