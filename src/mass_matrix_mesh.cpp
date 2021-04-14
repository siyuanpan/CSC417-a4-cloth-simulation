#include <mass_matrix_mesh.h>
#include <mass_matrix_linear_triangle.h>
#include <iostream>

void mass_matrix_mesh(Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::VectorXd> q,
                      Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F,
                      double density, Eigen::Ref<const Eigen::VectorXd> areas)
{
    typedef Eigen::Triplet<double> T;
    std::vector<T> triplets;
    triplets.reserve(F.rows() * 9 * 9);
    for (int i = 0; i < F.rows(); ++i)
    {
        Eigen::Matrix99d Mi;
        Eigen::RowVectorXi element = F.row(i);
        mass_matrix_linear_triangle(Mi, q, element, density, areas(i));
        for (int j = 0; j < 3; ++j)
        {
            for (int k = 0; k < 3; ++k)
            {
                for (int row = 0; row < 3; ++row)
                {
                    for (int col = 0; col < 3; ++col)
                    {
                        triplets.push_back(T(3 * element(j) + row, 3 * element(k) + col, Mi(3 * j + row, 3 * k + col)));
                    }
                }
            }
        }
    }
    M.resize(q.rows(), q.rows());
    M.setFromTriplets(triplets.begin(), triplets.end());
}
