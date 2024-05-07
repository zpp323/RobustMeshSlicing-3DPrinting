#include <main.h>

void checkFaceRobust(MyMesh& mesh, FaceHandle fh, float height) {
    float epsilon = 0.006;
    // cout << "height: " << height << endl;
    vector<MyMesh::Point> z;
    vector<VertexHandle> vhVector;
    for (MyMesh::FaceVertexIter fv_it = mesh.fv_iter(fh); fv_it.is_valid(); ++fv_it) {
        MyMesh::Point poi = mesh.point(*fv_it);
        z.push_back(poi);
        vhVector.push_back(*fv_it);
    }
    // cout << "z[0][2]: " << z[0][2] << endl;
    // cout << "z[1][2]: " << z[1][2] << endl;
    // cout << "z[2][2]: " << z[2][2] << endl;
    if (z[0][2] == height && z[0][2] == z[1][2] && z[1][2] == z[2][2]) {
        // cout << "condition 1" << endl;
        z[0][2] += epsilon;
        z[1][2] -= epsilon;
        mesh.set_point(vhVector[0], z[0]);
        mesh.set_point(vhVector[1], z[1]);
    }
    else if (z[0][2] == z[1][2] && z[0][2] == height) {
        // cout << "condition 2" << endl;
        z[0][2] += epsilon;
        z[1][2] -= epsilon;
        mesh.set_point(vhVector[0], z[0]);
        mesh.set_point(vhVector[1], z[1]);
    } else if (z[1][2] == z[2][2] && z[1][2] == height) {
        // cout << "condition 3" << endl;
        z[1][2] += epsilon;
        z[2][2] -= epsilon;
        mesh.set_point(vhVector[1], z[1]);
        mesh.set_point(vhVector[2], z[2]);
    } else if (z[2][2] == z[0][2] && z[2][2] == height) {
        // cout << "condition 4" << endl;
        z[2][2] += epsilon;
        z[0][2] -= epsilon;
        mesh.set_point(vhVector[2], z[2]);
        mesh.set_point(vhVector[0], z[0]);
    } else if (z[0][2] == height && z[1][2] > height && z[2][2] > height) {
        // cout << "condition 5" << endl;
        z[0][2] -= epsilon;
        mesh.set_point(vhVector[0], z[0]);
    } else if (z[0][2] == height && z[1][2] < height && z[2][2] < height) {
        // cout << "condition 6" << endl;
        z[0][2] += epsilon;
        mesh.set_point(vhVector[0], z[0]);
    } else if (z[1][2] == height && z[0][2] > height && z[2][2] > height) {
        // cout << "condition 7" << endl;
        z[1][2] -= epsilon;
        mesh.set_point(vhVector[1], z[1]);
    } else if (z[1][2] == height && z[0][2] < height && z[2][2] < height) {
        // cout << "condition 8" << endl;
        z[1][2] += epsilon;
        mesh.set_point(vhVector[1], z[1]);
    } else if (z[2][2] == height && z[0][2] > height && z[1][2] > height) {
        // cout << "condition 9" << endl;
        z[2][2] -= epsilon;
        mesh.set_point(vhVector[2], z[2]);
    } else if (z[2][2] == height && z[0][2] < height && z[1][2] < height) {
        // cout << "condition 10" << endl;
        z[2][2] += epsilon;
        mesh.set_point(vhVector[2], z[2]);
    }
}

void meshSlicing(MyMesh& mesh, float layerLevel, OpenMesh::FaceHandle fh, std::vector<float>& t,
std::vector<OpenMesh::FaceHandle>& sliceFaces, std::vector<OpenMesh::HalfedgeHandle>& sliceEdges) {
    int originFaceIdx = fh.idx();
    checkFaceRobust(mesh, fh, layerLevel);
    // cout << "check over" << endl;
    OpenMesh::FaceHandle next = fh;
    OpenMesh::FaceHandle nextNext;
    sliceFaces.push_back(fh);
    for (MyMesh::FaceHalfedgeIter fh_it = mesh.fh_iter(fh); fh_it.is_valid(); ++fh_it) {
        MyMesh::HalfedgeHandle heh = fh_it.handle();
        MyMesh::VertexHandle vh1 = mesh.from_vertex_handle(heh);
        MyMesh::VertexHandle vh2 = mesh.to_vertex_handle(heh);
        float z1 = mesh.point(vh1)[2];
        float z2 = mesh.point(vh2)[2];
        if (z1 >= layerLevel && z2 < layerLevel) {
            sliceEdges.push_back(heh);
            float t0 = (z1-layerLevel)/(z1-z2);
            t.push_back(t0);
            MyMesh::HalfedgeHandle heh2 = mesh.next_halfedge_handle(heh);
            float z3 = mesh.point(mesh.from_vertex_handle(heh2))[2];
            float z4 = mesh.point(mesh.to_vertex_handle(heh2))[2];
            if (z3 < layerLevel && z4 >= layerLevel) {
                t0 = (layerLevel-z3)/(z4-z3);
                t.push_back(t0);
                MyMesh::HalfedgeHandle opposite_heh = mesh.opposite_halfedge_handle(heh2);
                if (mesh.is_valid_handle(opposite_heh)) {  
                    MyMesh::FaceHandle fh_opposite = mesh.face_handle(opposite_heh);
                    next = fh_opposite;
                } else {
                    cout << "opposite edge handle error" << endl;
                }
            } else {
                heh2 = mesh.next_halfedge_handle(heh2);
                z3 = mesh.point(mesh.from_vertex_handle(heh2))[2];
                z4 = mesh.point(mesh.to_vertex_handle(heh2))[2];
                t0 = (layerLevel-z3)/(z4-z3);
                t.push_back(t0);
                MyMesh::HalfedgeHandle opposite_heh = mesh.opposite_halfedge_handle(heh2);
                if (mesh.is_valid_handle(opposite_heh)) {
                    MyMesh::FaceHandle fh_opposite = mesh.face_handle(opposite_heh);
                    next = fh_opposite;
                } else {
                    cout << "opposite edge handle error" << endl;
                }
            }
            sliceEdges.push_back(heh2);
            break;
        }
    }
    // cout << "this face: " << fh.idx() << endl;
    nextNext = next;
    while (next != fh) {
        checkFaceRobust(mesh, next, layerLevel);
        sliceFaces.push_back(next);
        // cout << "next face: " << next.idx() << endl;
        for (MyMesh::FaceHalfedgeIter fh_it = mesh.fh_iter(next); fh_it.is_valid(); ++fh_it) {
            MyMesh::HalfedgeHandle heh = fh_it.handle();
            MyMesh::VertexHandle vh1 = mesh.from_vertex_handle(heh);
            MyMesh::VertexHandle vh2 = mesh.to_vertex_handle(heh);
            float z1 = mesh.point(vh1)[2];
            float z2 = mesh.point(vh2)[2];
            // cout << "z1: " << z1 << endl;
            // cout << "z2: " << z2 << endl;
            if (z1 < layerLevel && z2 >= layerLevel) {
                float t1 = (layerLevel-z1)/(z2-z1);
                t.push_back(t1);
                MyMesh::HalfedgeHandle opposite_heh = mesh.opposite_halfedge_handle(heh);
                if (mesh.is_valid_handle(opposite_heh)) {
                    MyMesh::FaceHandle fh_opposite = mesh.face_handle(opposite_heh);
                    next = fh_opposite;
                } else {
                    cout << "opposite edge handle error" << endl;
                }
                sliceEdges.push_back(heh);
                break;
            }
        }
        if (next == nextNext) {
            cout << "error" << endl;
            break;
        }
    }
    
    
}

int main() {
    bool is_tet_wild = false;
    string inputF;
    cout << "Please input the file name:" << endl;
    cin >> inputF;
    const char* inputFileName = inputF.c_str();
    char sliceFileName[100] = "mid.obj";
    char tempFileName[100] = {0};
    strcpy(tempFileName, inputFileName);
    pcl::PointCloud<pcl::PointXYZ> ::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
    char* pext = strrchr(tempFileName, '.');
    std::string extply("ply");
    std::string extpcd("pcd");
    if (pext){
        *pext = '\0';
        pext++;
    }
    std::string ext(pext);
    if (!((ext == extply) || (ext == extpcd))){
        is_tet_wild = true;
    } else {    // poisson surface reconstruction
        cout << "execute poisson surface reconstruction" << endl;
        if (ext == extply){
            if (pcl::io::loadPLYFile(inputFileName, *cloud) == -1){
                PCL_ERROR("Could not read ply file!\n");
                return -1;
            }
        }
        else{
            if (pcl::io::loadPCDFile(inputFileName, *cloud) == -1){
                PCL_ERROR("Could not read pcd file!\n");
                return -1;
            }
        }
        cout << "point could size: " << cloud->points.size() << endl;

        pcl::PointCloud<pcl::PointNormal>::Ptr cloud_with_normals(new pcl::PointCloud<pcl::PointNormal>); //法向量点云对象指针
        pcl::NormalEstimation<pcl::PointXYZ , pcl::Normal> n;
        pcl::PointCloud<pcl::Normal>::Ptr normals(new pcl::PointCloud<pcl::Normal>);
        pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>);
        tree->setInputCloud(cloud);
        n.setInputCloud(cloud);
        n.setSearchMethod(tree);
        n.setKSearch(20);
        n.compute(*normals);

        pcl::concatenateFields(*cloud , *normals , *cloud_with_normals);

        pcl::search::KdTree<pcl::PointNormal>::Ptr tree2(new pcl::search::KdTree<pcl::PointNormal>) ;
        tree2->setInputCloud(cloud_with_normals);

        pcl::Poisson<pcl::PointNormal> pn;
        pn.setConfidence(false);
        pn.setDegree(2);
        pn.setDepth(8);
        pn.setIsoDivide(8);
        pn.setManifold(false);
        pn.setOutputPolygons(false);
        pn.setSamplesPerNode(9);
        pn.setScale(1.25);
        pn.setSolverDivide(8);

        pn.setSearchMethod(tree2);
        pn.setInputCloud(cloud_with_normals);
        pcl::PolygonMesh p_mesh;
        //执行重构
        pn.performReconstruction(p_mesh);

        //保存网格图
        pcl::io::saveOBJFile(sliceFileName, p_mesh);
        cout << "construction ok" << endl;
    }
    // tetwild
    if (is_tet_wild) {
        cout << "execute poisson surface reconstruction" << endl;
        Eigen::MatrixXd v_in;
        Eigen::MatrixXi f_in;
        Eigen::MatrixXd v_out;
        Eigen::MatrixXi t_out;
        Eigen::VectorXd a_out;

        Eigen::MatrixXd v_s;
        Eigen::MatrixXi f_s;
        igl::readOBJ(inputFileName, v_in, f_in);
        tetwild::Args args;
        tetwild::tetrahedralization(v_in, f_in, v_out, t_out, a_out, args);
        tetwild::extractSurfaceMesh(v_out, t_out, v_s, f_s);
        igl::writeOBJ(sliceFileName, v_s, f_s);
    }

    /* slicing a mesh */

    // read mesh
    MyMesh mesh;
    if (!OpenMesh::IO::read_mesh(mesh, sliceFileName)) {
        std::cerr << "read error\n";
        exit(1);
    }
    cout << "read mesh ok!" << endl;
    
    // mesh slicing
    float zMin = INF;
    float zMax = (-1)*INF;

    // sorting
    std::vector<faceCompare> faces;
    std::vector<float> t;
    std::vector<OpenMesh::FaceHandle> sliceFaces;
    std::vector<OpenMesh::HalfedgeHandle> sliceEdges;
    for (OpenMesh::TriMesh_ArrayKernelT<>::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it) {
        float min_z = INF;
        float max_z = (-1)*INF;
        for (MyMesh::FaceVertexIter it = mesh.fv_begin(*f_it); it.is_valid(); ++it) {
	        MyMesh::Point myPoint = mesh.point(*it);
            if (min_z > myPoint[2]) {
                min_z = myPoint[2];
            }
            if (max_z < myPoint[2]) {
                max_z = myPoint[2];
            }
        }
        faceCompare fc;
        fc.fh = f_it.handle();
        fc.min_z = min_z;
        fc.max_z = max_z;
        faces.push_back(fc);
        zMin = std::min(zMin, min_z);
        zMax = std::max(zMax, max_z);
    }
    std::sort(faces.begin(), faces.end(), [](faceCompare& lhs, faceCompare& rhs){
            if (lhs.min_z != rhs.min_z) {
                return lhs.min_z < rhs.min_z;
            }
            return lhs.max_z < rhs.max_z;
        });
    // for (auto ele: faces) {
    //     cout << ele.min_z << " " << ele.max_z << endl;
    // }

    // slicing
    int cursor = 0;
    float layerLevel, layerBase;
    int layerNumber = 0;
    float layerThick = 1;
    layerLevel = layerBase = zMin + layerThick / 10;
    while (layerLevel < zMax - layerThick / 10) {
        layerNumber++;
        cout << "layer height: " << layerLevel << endl;
        t.clear();
        sliceFaces.clear();
        sliceEdges.clear();
        while (cursor<faces.size() && faces[cursor].min_z<layerLevel && faces[cursor].max_z<layerLevel) {
            cursor++;
        }
        if (cursor == faces.size()) {
            cout << "over" << endl;
            break;
        }
        // cout << faces.size() << " " << cursor << " min height: " << faces[cursor].min_z << " max height: " << faces[cursor].max_z << " layer height: " << layerLevel << endl;
        // triangle-plane intersection
        meshSlicing(mesh, layerLevel, faces[cursor].fh, t, sliceFaces, sliceEdges);
        t.pop_back();
        // cout << "t: " << endl;
        // for (int i=0;i<t.size();i++){
        //     cout << t[i] << endl;
        // }
        // cout << sliceFaces.size() << endl;
        layerLevel += layerThick;
        // output the contour
        sliceEdges.pop_back();
        ofstream ofs;
        ofs.open("data/contour_"+to_string(layerNumber)+".obj");
        for (int i=0;i<t.size();i++) {
            MyMesh::VertexHandle vh1 = mesh.from_vertex_handle(sliceEdges[i]);
            MyMesh::VertexHandle vh2 = mesh.to_vertex_handle(sliceEdges[i]);
            MyMesh::Point p1 = mesh.point(vh1);
            MyMesh::Point p2 = mesh.point(vh2);
            MyMesh::Point p3 = p2-p1;
            MyMesh::Point pMid = p1 + t[i] * p3;
            // pMid[0] = p1[0] + t[i] * p3[0];
            // pMid[1] = p1[1] + t * p3[1];
            // pMid[2] = p1[2] + t * p3[2];
            // cout << "x: " << pMid[0] << " y: " << pMid[1] << " z: " << pMid[2] << endl;
            ofs << "v " << pMid[0] << " " << pMid[1] << " " << pMid[2] << endl;
        }
        for (int i=0;i<t.size()-1;i++) {
            ofs << "l " << i+1 << " " << i+2 << endl;
        }
        ofs << "l " << t.size() << " " << 1 << endl;
        ofs.close();
    }
    cout << "slicing ok! " << layerNumber << " layers in total." << endl;
    return 0;
}