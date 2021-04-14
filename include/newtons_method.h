#include <Eigen/Dense>
#include <EigenTypes.h>

//Input:
//  x0 - initial point for newtons search
//  f(x) - function that evaluates and returns the cost function at x
//  g(dx, x) - function that evaluates and returns the gradient of the cost function in dx
//  H(dH, x) - function that evaluates and returns the Hessian in dH (as a sparse matrix).
//  max steps - the maximum newton iterations to take
//  tmp_g and tmp_H are scratch space to store gradients and hessians
//Output:
//  x0 - update x0 to new value
template <typename Objective, typename Jacobian, typename Hessian>
double newtons_method(Eigen::VectorXd &x0, Objective &f, Jacobian &g, Hessian &H, unsigned int maxSteps, Eigen::VectorXd &tmp_g, Eigen::SparseMatrixd &tmp_H)
{
    for (unsigned int i = 0; i < maxSteps; ++i)
    {
        g(tmp_g, x0);
        // std::cout << tmp_g.squaredNorm() << std::endl;
        if (tmp_g.squaredNorm() < std::numeric_limits<double>::epsilon())
            return 0.0;
        Eigen::SimplicialLDLT<Eigen::SparseMatrixd, Eigen::Upper> solver;
        // Eigen::SimplicialLLT<Eigen::SparseMatrixd> solver;
        H(tmp_H, x0);
        solver.analyzePattern(tmp_H);
        solver.factorize(tmp_H);

        double Regularization = 0.00001;
        bool success = true;
        Eigen::SparseMatrix<double> I(tmp_H.rows(), tmp_H.rows());

        while (solver.info() != Eigen::Success)
        {
            Regularization *= 10;
            tmp_H = tmp_H + Regularization * I;
            solver.factorize(tmp_H);
            success = false;
            // std::cout << "solver fail\n";
            // return 0.0;
        }
        if (!success)
            std::cout << " adding " << Regularization
                      << " identites.(ldlt solver)" << std::endl;

        // std::cout << tmp_g << std::endl;
        // std::exit(1);
        Eigen::VectorXd d;
        d = -solver.solve(tmp_g);
        // std::cout << d(0) << std::endl;
        // std::exit(0);

        double alpha = 1;
        double p = 0.5;
        double c = 1e-8;

        for (;;)
        {
            if (f(x0 + alpha * d) <= f(x0) + c * d.transpose() * tmp_g)
                break;
            alpha *= p;
            if (alpha < std::numeric_limits<double>::epsilon())
                return 0.0;
        }

        x0 += alpha * d;
    }

    return 0.0;
}
