cmake_minimum_required(VERSION 3.0)

project(fastSlicing LANGUAGES CXX)

include_directories(.)
include_directories(../TetWild/include)
include_directories(../TetWild/extern/libigl/include)

link_directories (../TetWild/build)
link_directories (../TetWild/build/lib)
link_directories (../TetWild/build/extern/fmt)
link_directories (../TetWild/build/extern/pymesh)

set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -O3  -Wall")

find_package(Boost REQUIRED
             COMPONENTS system filesystem thread)

find_package(PCL REQUIRED)
include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

add_executable(fastSlicing
        main.cpp main.h cloud.h)

target_link_libraries (${PROJECT_NAME} OpenMeshCore OpenMeshTools)

target_link_libraries (${PROJECT_NAME} tetwild fmt pymesh_tiny gmp mpfr geogram boost_thread)

target_link_libraries (${PROJECT_NAME}  ${PCL_LIBRARIES})
