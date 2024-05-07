#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/Handles.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/Traits.hh>
#include <OpenMesh/Core/System/config.h>
#include <OpenMesh/Core/Mesh/DefaultTriMesh.hh>
#include <OpenMesh/Core/Utils/Property.hh>
#include <OpenMesh/Core/Utils/PropertyManager.hh>

#include <tetwild/tetwild.h>

#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>

#include "cloud.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>

using namespace std;
using namespace OpenMesh;

#define INF 1e8  

typedef OpenMesh::PolyMesh_ArrayKernelT<>  MyMesh;  

struct faceCompare {
    OpenMesh::FaceHandle fh;
    float min_z;
    float max_z;
};
