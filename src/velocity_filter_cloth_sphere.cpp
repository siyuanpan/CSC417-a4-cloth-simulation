#include <velocity_filter_cloth_sphere.h>

void velocity_filter_cloth_sphere(Eigen::VectorXd &qdot, const std::vector<unsigned int> &indices, 
                                  const std::vector<Eigen::Vector3d> &normals) {
    for (int i = 0; i < indices.size(); ++i) {
        auto idx = indices[i];
        Eigen::Vector3d v = qdot.segment<3>(3*idx);
        Eigen::Vector3d n = normals[i];
        double a;
        if (n.transpose()*v >= 0) {
            a = 0;
        } else {
            a = -n.transpose()*v;
        }
        qdot.segment<3>(3*idx) = v + a*n;
    }

    

}
