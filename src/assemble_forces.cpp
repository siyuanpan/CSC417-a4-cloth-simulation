#include <assemble_forces.h>
#include <dV_membrane_corotational_dq.h>
#include <dphi_cloth_triangle_dX.h>
#include <iostream>

void assemble_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::MatrixXd> qdot, Eigen::Ref<const Eigen::MatrixXd> dX,
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F, Eigen::Ref<const Eigen::VectorXd> a0,
                     double mu, double lambda) {
    f.resize(q.rows());
    f.setZero();

#define force(i, j) f.segment<3>(3 * F((i), (j)))
    for (int i = 0; i < F.rows(); ++i)
    {
        Eigen::Vector9d dV;
        dV_membrane_corotational_dq(dV, q, dX, V, F.row(i), a0(i), mu, lambda);
        force(i, 0) += -dV.segment<3>(0);
        force(i, 1) += -dV.segment<3>(3);
        force(i, 2) += -dV.segment<3>(6);
    }
#undef force
        
        

};
