#include <assemble_forces.h>
#include <dV_membrane_corotational_dq.h>
#include <dphi_cloth_triangle_dX.h>
#include <iostream>
#include <dV_spring_particle_particle_dq.h>
#include <igl/edges.h>
#include <igl/edge_lengths.h>

void assemble_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::MatrixXd> qdot, Eigen::Ref<const Eigen::MatrixXd> dX,
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F, Eigen::Ref<const Eigen::VectorXd> a0,
                     double mu, double lambda)
{
    f.resize(q.rows());
    f.setZero();

#define force(i, j) f.segment<3>(3 * F((i), (j)))
    // Eigen::MatrixXi E;
    // Eigen::VectorXd L;
    // igl::edges(F, E);
    // igl::edge_lengths(V, E, L);

    // for (int i = 0; i < E.rows(); ++i)
    // {
    //     auto q0 = q.segment<3>(3 * E(i, 0));
    //     auto q1 = q.segment<3>(3 * E(i, 1));
    //     Eigen::Vector6d fi;
    //     dV_spring_particle_particle_dq(fi, q0, q1, L(i), 100);
    //     f.segment<3>(3 * E(i, 0)) -= fi.segment<3>(0);
    //     f.segment<3>(3 * E(i, 1)) -= fi.segment<3>(3);
    // }

    for (int i = 0; i < F.rows(); ++i)
    {
        Eigen::Vector9d dV;
        Eigen::Vector9d tmp = dX.row(i);
        Eigen::Matrix3d dx = Eigen::Map<Eigen::MatrixXd>(tmp.data(), 3, 3);
        // Eigen::Map<Eigen::Matrix3d>(dX.row(i).data())
        dV_membrane_corotational_dq(dV, q, dx, V, F.row(i), a0(i), mu, lambda);
        // std::cout << dV << std::endl;
        force(i, 0) += -dV.segment<3>(0);
        force(i, 1) += -dV.segment<3>(3);
        force(i, 2) += -dV.segment<3>(6);
    }
#undef force

    //std::cout << f.transpose() << std::endl;
    // std::exit(1);
};
