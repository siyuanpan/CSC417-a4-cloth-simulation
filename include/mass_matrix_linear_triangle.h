#include <Eigen/Dense>
#include <EigenTypes.h>


void mass_matrix_linear_triangle(Eigen::Matrix99d &M, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::RowVectorXi> element, double density, double volume);
