# RobustMeshSlicing-3DPrinting



V1.0

Input file format: obj/ply/pcd

requirement:

1. [tetWild](https://github.com/Yixin-Hu/TetWild) : This library needs to be compiled and the include directories and link directories in CMakeLists.txt in this project needs to change with it.

```cmake
include_directories(../TetWild/include)
include_directories(../TetWild/extern/libigl/include)

link_directories (../TetWild/build)
link_directories (../TetWild/build/lib)
link_directories (../TetWild/build/extern/fmt)
link_directories (../TetWild/build/extern/pymesh)
```

2. [OpenMesh](https://www.graphics.rwth-aachen.de/software/openmesh/)
3. [pcl](https://pointclouds.org): This library could be installed in this way in Ubuntu:

```shell
sudo apt-get install libpcl-dev
```