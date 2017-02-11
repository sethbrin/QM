cmake_minimum_required(VERSION 2.8)
project(version2)

set(CXX_FLAGS
    -g
    -std=c++11
    -O0
    )

string(REPLACE ";" " " CMAKE_CXX_FLAGS "${CXX_FLAGS}")

set(LIB_SOURCE_FILES
    gzstream/gzstream.h
    common.cpp
    common.h
    grid_structures.cpp
    grid_structures.h
    qm_interpolation.cpp
    qm_interpolation.h
    read_energy_force_new.cpp
    read_energy_force_new.h
    )

set(SOURCE_FILES
    ${LIB_SOURCE_FILES}
    main.cpp)


link_directories(${PROJECT_SOURCE_DIR}/gzstream)

add_library(waterbox ${LIB_SOURCE_FILES})
add_executable(version2 ${SOURCE_FILES})
target_link_libraries(version2
    z gzstream.a boost_filesystem boost_system boost_thread)

