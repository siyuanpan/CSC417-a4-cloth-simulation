#include <dV_spring_particle_particle_dq.h>

void dV_spring_particle_particle_dq(Eigen::Ref<Eigen::Vector6d> f, Eigen::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d>     q1, double l0, double stiffness) {
    auto k = stiffness*((q0-q1).norm()-l0);
    auto direction = (q0-q1).normalized();
    f.segment<3>(0) = k * direction;
    f.segment<3>(3) = -k * direction;
   
    
}
