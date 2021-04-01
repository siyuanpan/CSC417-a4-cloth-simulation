#include <Eigen/Dense>
#include <Eigen/Sparse>
#include<Eigen/SparseCholesky>
#include <EigenTypes.h>

//Input:
//  q - generalized coordinates for the FEM system
//  qdot - generalized velocity for the FEM system
//  dt - the time step in seconds
//  mass - the mass matrix
//  force(f, q, qdot) - a function that computes the force acting on the FEM system. This takes q and qdot as parameters, returns the force in f.
//  stiffness(K, q, qdot) - a function that computes the stiffness (negative second derivative of the potential energy). This takes q and qdot as parameters, returns the stiffness matrix in K.  
//  tmp_force - scratch space to collect forces
//  tmp_stiffness - scratch space to collect stiffness matrix
//Output:
//  q - set q to the updated generalized coordinate using linearly implicit time integration
//  qdot - set qdot to the updated generalized velocity using linearly implicit time integration
template<typename FORCE, typename STIFFNESS> 
inline void linearly_implicit_euler(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, 
                            const Eigen::SparseMatrixd &mass,  FORCE &force, STIFFNESS &stiffness, 
                            Eigen::VectorXd &tmp_force, Eigen::SparseMatrixd &tmp_stiffness) {
     stiffness(tmp_stiffness, q, qdot);
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    tmp_stiffness = (mass - dt * dt * tmp_stiffness).eval();
    solver.analyzePattern(tmp_stiffness);
    solver.factorize(tmp_stiffness);
    // std::cout << tmp_stiffness.block(0, 0, 32, 32) << std::endl;
    // std::exit(0);
    double Regularization = 0.00001;
    bool success = true;
    Eigen::SparseMatrix<double> I(tmp_stiffness.rows(), tmp_stiffness.rows());

    if (solver.info() != Eigen::Success)
    {
        Regularization *= 10;
        tmp_stiffness = tmp_stiffness + Regularization * I;
        solver.factorize(tmp_stiffness);
        success = false;
    }
    if (!success)
        std::cout << " adding " << Regularization
                  << " identites.(ldlt solver)" << std::endl;

    force(tmp_force, q, qdot);
    tmp_force = (mass * qdot + dt * tmp_force).eval();

    qdot = solver.solve(tmp_force);
    q = (q + dt * qdot).eval();
    


}
