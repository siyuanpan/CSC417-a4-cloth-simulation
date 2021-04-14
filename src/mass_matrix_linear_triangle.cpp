#include "mass_matrix_linear_triangle.h"

void mass_matrix_linear_triangle(Eigen::Matrix99d &M, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::RowVectorXi> element, double density, double volume)
{
    double mass = density * volume;
    double c0 = (1.0 / 6.0) * mass;
    double c1 = (1.0 / 12.0) * mass;
    // double c0 = (1.0 / 10.0) * mass;
    // double c1 = (1.0 / 20.0) * mass;

    M.block(0, 0, 3, 3) = c0 * Eigen::Matrix3d::Identity();
    M.block(0, 3, 3, 3) = c1 * Eigen::Matrix3d::Identity();
    M.block(0, 6, 3, 3) = c1 * Eigen::Matrix3d::Identity();

    M.block(3, 0, 3, 3) = c1 * Eigen::Matrix3d::Identity();
    M.block(3, 3, 3, 3) = c0 * Eigen::Matrix3d::Identity();
    M.block(3, 6, 3, 3) = c1 * Eigen::Matrix3d::Identity();

    M.block(6, 0, 3, 3) = c1 * Eigen::Matrix3d::Identity();
    M.block(6, 3, 3, 3) = c1 * Eigen::Matrix3d::Identity();
    M.block(6, 6, 3, 3) = c0 * Eigen::Matrix3d::Identity();
}
