// #include "stdafx.h"
#include <iostream>
#include <string>

#include <pcl/io/pcd_io.h>
#include <pcl/io/ply_io.h>
#include <pcl/io/obj_io.h>

#include <pcl/point_types.h>

#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/features/normal_3d.h>
#include <pcl/features/normal_3d_omp.h>

#include <pcl/surface/poisson.h>
#include <pcl/surface/gp3.h>

#include <pcl/visualization/cloud_viewer.h>
#include <pcl/visualization/pcl_visualizer.h>

//多线程
#include <boost/thread/thread.hpp>

#include <vector>
using namespace std;

struct Node
{
	float x;
	float y;
	float z;
	float i;
	float j;
	float k;
};
