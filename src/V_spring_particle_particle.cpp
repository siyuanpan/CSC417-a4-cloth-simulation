#include <V_spring_particle_particle.h>

void V_spring_particle_particle(double &V, Eigen ::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d> q1, double l0, double stiffness) {
    V = 0.5*stiffness*std::pow((q0-q1).norm()-l0, 2);
    
    
}
