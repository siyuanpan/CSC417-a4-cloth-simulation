#include <dV_membrane_corotational_dq.h>
#include <dphi_cloth_triangle_dX.h>
#include <iostream>

void dV_membrane_corotational_dq(Eigen::Vector9d &dV, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::Matrix3d> dX,
                                 Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double area,
                                 double mu, double lambda)
{

    //Deformation Gradient
    Eigen::Matrix3d dx; //deformed tangent matrix
    Eigen::Matrix3d U;
    Eigen::Vector3d S;
    Eigen::Matrix3d W;

    //TODO: SVD Here
    Eigen::Vector3d x0 = q.segment<3>(3 * element(0));
    Eigen::Vector3d x1 = q.segment<3>(3 * element(1));
    Eigen::Vector3d x2 = q.segment<3>(3 * element(2));
    Eigen::Vector3d n = (x1 - x0).cross(x2 - x0);

    Eigen::Matrix<double, 3, 4> left;
    left.col(0) = x0;
    left.col(1) = x1;
    left.col(2) = x2;
    left.col(3) = n.normalized();

    Eigen::Vector3d X0 = V.row(element(0)).transpose();
    Eigen::Vector3d X1 = V.row(element(1)).transpose();
    Eigen::Vector3d X2 = V.row(element(2)).transpose();
    Eigen::Vector3d N = (X1 - X0).cross(X2 - X0).normalized();
    Eigen::Matrix<double, 4, 3> right;
    // dphi_cloth_triangle_dX(dx, V, element, X0);
    right.block<3, 3>(0, 0) = dX;
    right.row(3) = N.transpose();
    Eigen::Matrix3d F = left * right;
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    U = svd.matrixU();
    S = svd.singularValues();
    W = svd.matrixV();
    // std::cout << F << std::endl;
    // std::exit(1);

    //Fix for inverted elements (thanks to Danny Kaufman)
    double det = S[0] * S[1];

    if (det <= -1e-10)
    {
        if (S[0] < 0)
            S[0] *= -1;
        if (S[1] < 0)
            S[1] *= -1;
        if (S[2] < 0)
            S[2] *= -1;
    }

    if (U.determinant() <= 0)
    {
        U(0, 2) *= -1;
        U(1, 2) *= -1;
        U(2, 2) *= -1;
    }

    if (W.determinant() <= 0)
    {
        W(0, 2) *= -1;
        W(1, 2) *= -1;
        W(2, 2) *= -1;
    }

    // F = U * S.asDiagonal() * W.transpose();

    //TODO: energy model gradient
    Eigen::Matrix3d dS;
    dS << 2 * mu * (S[0] - 1) + lambda * (S[0] + S[1] + S[2] - 3), 0, 0,
        0, 2 * mu * (S[1] - 1) + lambda * (S[0] + S[1] + S[2] - 3), 0,
        0, 0, 2 * mu * (S[2] - 1) + lambda * (S[0] + S[1] + S[2] - 3);
    Eigen::Matrix3d dpsi_dF = U * dS * W.transpose();
    // std::cout << dS << std::endl;

    Eigen::Matrix<double, 9, 9> B;
    B.setZero();
    for (int i = 0; i < 3; ++i)
    {
        B.block(0, 0 + 3 * i, 3, 1) = dX.row(i).transpose();
        B.block(3, 1 + 3 * i, 3, 1) = dX.row(i).transpose();
        B.block(6, 2 + 3 * i, 3, 1) = dX.row(i).transpose();
    }
    Eigen::Matrix<double, 9, 3> N_;
    N_.setZero();
    N_.block(0, 0, 3, 1) = N;
    N_.block(3, 1, 3, 1) = N;
    N_.block(6, 2, 3, 1) = N;

    auto crossMat = [](Eigen::Matrix3d &mat, Eigen::Vector3d &v1, Eigen::Vector3d &v2) {
        Eigen::Vector3d v = v2 - v1;
        mat << 0, -v[2], v[1],
            v[2], 0, -v[0],
            -v[1], v[0], 0;
    };
    Eigen::Matrix3d dx1, dx2;
    crossMat(dx1, x0, x1);
    crossMat(dx2, x0, x2);
    Eigen::Matrix<double, 3, 9> tmp1, tmp2;
    tmp1.setZero();
    tmp1.block(0, 0, 3, 3) = -Eigen::Matrix3d::Identity();
    tmp1.block(0, 6, 3, 3) = Eigen::Matrix3d::Identity();
    tmp2.setZero();
    tmp2.block(0, 0, 3, 3) = -Eigen::Matrix3d::Identity();
    tmp2.block(0, 3, 3, 3) = Eigen::Matrix3d::Identity();

    Eigen::Matrix<double, 9, 9> dF;
    Eigen::Vector3d nn = n.normalized();
    dF = B + N_ * (Eigen::Matrix3d::Identity() - nn * nn.transpose()) * (dx1 * tmp1 - dx2 * tmp2) / n.norm();

    Eigen::Matrix3d psiT = dpsi_dF.transpose();
    Eigen::Vector9d flat = Eigen::Map<Eigen::Vector9d>(psiT.data(), 9, 1);
    // std::cout << dF << std::endl;
    // std::cout << dpsi_dF << std::endl;
    // std::cout << psiT << std::endl;
    // std::cout << flat << std::endl;
    // std::exit(1);

    dV = area * dF.transpose() * flat;
    // std::cout << dV << std::endl;
    // std::exit(1);
}
