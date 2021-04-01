#include <assemble_stiffness.h>
#include <d2V_membrane_corotational_dq2.h>
#include <iostream>
void assemble_stiffness(Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::MatrixXd> dX,
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F, Eigen::Ref<const Eigen::VectorXd> a0, 
                     double mu, double lambda) {
    typedef Eigen::Triplet<double> T;
    std::vector<T> triplets;
    triplets.reserve(F.rows() * 4 * 4 * 9);
    for (int i = 0; i < F.rows(); ++i)
    {
        Eigen::Matrix99d dV;
        d2V_membrane_corotational_dq2(dV, q, dX, V, F.row(i), a0(i), mu, lambda);

        for (int j = 0; j < 3; ++j)
        {
            for (int k = 0; k < 3; ++k)
            {
                for (int row = 0; row < 3; ++row)
                {
                    for (int col = 0; col < 3; ++col)
                    {
                        triplets.push_back(T(3 * F(i, j) + row, 3 * F(i, k) + col, -dV(3 * j + row, 3 * k + col)));
                    }
                }
            }
        }
    }
    K.resize(q.rows(), q.rows());
    K.setFromTriplets(triplets.begin(), triplets.end());
        
       
        
};
