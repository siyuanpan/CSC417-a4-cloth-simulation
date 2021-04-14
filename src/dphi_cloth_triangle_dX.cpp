#include <dphi_cloth_triangle_dX.h>

//compute 3x3 deformation gradient
void dphi_cloth_triangle_dX(Eigen::Matrix3d &dphi, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X)
{
    auto x0 = V.row(element(0)).transpose();
    auto x1 = V.row(element(1)).transpose();
    auto x2 = V.row(element(2)).transpose();

    Eigen::Matrix32d T;
    T.col(0) = x1 - x0;
    T.col(1) = x2 - x0;

    Eigen::Matrix<double, 2, 3> TT = T.transpose();
    Eigen::Vector2d NegOne = -Eigen::Vector2d::Ones();
    dphi.row(0) = NegOne.transpose() * (TT * T).inverse() * TT;
    dphi.block<2, 3>(1, 0) = (TT * T).inverse() * TT;
    // dphi.noalias() = dphi.transpose();
}
