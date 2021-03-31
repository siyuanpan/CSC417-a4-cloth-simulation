#include <V_membrane_corotational.h>

//Allowed to use libigl SVD or Eigen SVD for this part
void V_membrane_corotational(double &energy, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::Matrix3d> dX, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double area, 
                          double mu, double lambda) {
    auto x0 = q.segment<3>(3*element(0));
    auto x1 = q.segment<3>(3*element(1));
    auto x2 = q.segment<3>(3*element(2));
    auto n = (x1-x0).cross(x2-x0).normalized();

    Eigen::Matrix<double, 3, 4> left;
    left.col(0) = x0;
    left.col(1) = x1;
    left.col(2) = x2;
    left.col(3) = n;

    Eigen::Vector3d X0 = V.row(element(0)).transpose();
    Eigen::Vector3d X1 = V.row(element(1)).transpose();
    Eigen::Vector3d X2 = V.row(element(2)).transpose();
    // Eigen::Vector3d N = (V.row(element(1)).transpose() - V.row(element(0)).transpose()).cross(V.row(element(2)).transpose() - V.row(element(0)).transpose());
    Eigen::Vector3d N = (X1-X0).cross(X2-X0).normalized();
    Eigen::Matrix<double, 4, 3> right;
    right.block<3, 3>(0, 0) = dX;
    right.row(3) = N.transpose();

    Eigen::Matrix3d F = left * right;
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(F);

    auto sigma = svd.singularValues();
    energy = area * mu * (std::pow(sigma(0)-1, 2) + std::pow(sigma(1)-1, 2) + std::pow(sigma(2)-1, 2)) + 0.5*lambda*std::pow(sigma(0)+sigma(1)+sigma(2)-3, 2);
}
