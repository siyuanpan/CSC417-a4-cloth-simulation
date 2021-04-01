#include <d2V_membrane_corotational_dq2.h>
#include <dsvd.h>
#include <iostream>

void d2V_membrane_corotational_dq2(Eigen::Matrix99d &H, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::Matrix3d> dX, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double area, 
                          double mu, double lambda) {
    

    //SVD = USW^T
    Eigen::Matrix3d U;
    Eigen::Vector3d S; 
    Eigen::Matrix3d W; 
    Eigen::Matrix3d F; //deformation gradient
    
    double tol = 1e-5;
    
    //Compute SVD of F here
    Eigen::Vector3d x0 = q.segment<3>(3*element(0));
    Eigen::Vector3d x1 = q.segment<3>(3*element(1));
    Eigen::Vector3d x2 = q.segment<3>(3*element(2));
    Eigen::Vector3d n = (x1-x0).cross(x2-x0);

    Eigen::Matrix<double, 3, 4> left;
    left.col(0) = x0;
    left.col(1) = x1;
    left.col(2) = x2;
    left.col(3) = n.normalized();

    Eigen::Vector3d X0 = V.row(element(0)).transpose();
    Eigen::Vector3d X1 = V.row(element(1)).transpose();
    Eigen::Vector3d X2 = V.row(element(2)).transpose();
    Eigen::Vector3d N = (X1-X0).cross(X2-X0).normalized();
    Eigen::Matrix<double, 4, 3> right;
    //dphi_cloth_triangle_dX(dx, V, element, X0);
    right.block<3, 3>(0, 0) = dX;
    right.row(3) = N.transpose();
    F = left * right;
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    U = svd.matrixU();
    S = svd.singularValues();
    W = svd.matrixV();

    //deal with singularity in the svd gradient
    if(std::fabs(S[0] - S[1]) < tol || std::fabs(S[1] - S[2]) < tol || std::fabs(S[0] - S[2]) < tol) {
        F += Eigen::Matrix3d::Random()*tol;
        Eigen::JacobiSVD<Eigen::Matrix3d> svd2(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
        U = svd2.matrixU();
        W = svd2.matrixV();
        S = svd2.singularValues();
    }
    
    //Fix for inverted elements (thanks to Danny Kaufman)
    double det = S[0]*S[1];
    
     if(det <= -1e-10)
    {
        if(S[0] < 0) S[0] *= -1;
        if(S[1] < 0) S[1] *= -1;
        if(S[2] < 0) S[2] *= -1;
    }
    
    if(U.determinant() <= 0)
    {
        U(0, 2) *= -1;
        U(1, 2) *= -1;
        U(2, 2) *= -1;
    }
    
    if(W.determinant() <= 0)
    {
        W(0, 2) *= -1;
        W(1, 2) *= -1;
        W(2, 2) *= -1;
    }

    //TODO: compute H, the hessian of the corotational energy
    Eigen::Tensor3333d dU, dV;
    Eigen::Tensor333d dS;
    dsvd(dU, dS, dV, F);
    double parm1 = 2*mu+lambda;
    double parm2 = lambda;
    Eigen::Matrix3d dpsids;
    dpsids << 2*mu*(S[0]-1)+lambda*(S[0]+S[1]+S[2]-3), 0, 0,
        0, 2*mu*(S[1]-1)+lambda*(S[0]+S[1]+S[2]-3), 0,
        0, 0, 2*mu*(S[2]-1)+lambda*(S[0]+S[1]+S[2]-3);
    Eigen::Matrix3d d2psi;
    d2psi << parm1, parm2, parm2,
        parm2, parm1, parm2,
        parm2, parm2, parm1;
    Eigen::Matrix<double, 9, 9> d2psidF;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            Eigen::Matrix3d dUdF;
            Eigen::Matrix3d dVdF;
            Eigen::Vector3d dSdF;
            dSdF << dS[0][i][j], dS[1][i][j], dS[2][i][j];
            for (int x = 0; x < 3; ++x) {
                for (int y = 0; y < 3; ++y) {
                    dUdF(x, y) = dU[x][y](i, j);
                    dVdF(x, y) = dV[x][y](i, j);
                }
            }
            Eigen::Vector3d dsij = d2psi*dSdF;
            Eigen::Matrix<double, 3, 3, Eigen::RowMajor> d2psidFij = dUdF*dpsids*W.transpose() + U*dsij.asDiagonal()*W.transpose() + U*dpsids*dVdF.transpose();
            Eigen::Vector9d tmp = Eigen::Map<Eigen::Vector9d>(d2psidFij.data(), 9);
            d2psidF.col(i*3+j) = tmp;
        }
    }

    Eigen::Matrix<double, 9, 9> B;
    B.setZero();
    for (int i = 0; i < 3; ++i) {
        B.block(0, 0+3*i, 3, 1) = dX.row(i).transpose();
        B.block(3, 1+3*i, 3, 1) = dX.row(i).transpose();
        B.block(6, 2+3*i, 3, 1) = dX.row(i).transpose();
    }
    Eigen::Matrix<double, 9, 3> N_;
    N_.setZero();
    N_.block(0, 0, 3, 1) = N;
    N_.block(3, 1, 3, 1) = N;
    N_.block(6, 2, 3, 1) = N;

    auto crossMat = [](Eigen::Matrix3d& mat, Eigen::Vector3d& v1, Eigen::Vector3d& v2) {
        Eigen::Vector3d v = v2-v1;
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
    dF = B + N_ * (Eigen::Matrix3d::Identity()-nn*nn.transpose()) * (dx1*tmp1-dx2*tmp2)/n.norm();

    H = dF.transpose()*d2psidF*dF;

    //fix errant eigenvalues
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix99d> es(H);
    
    Eigen::MatrixXd DiagEval = es.eigenvalues().real().asDiagonal();
    Eigen::MatrixXd Evec = es.eigenvectors().real();
    
    for (int i = 0; i < 9; ++i) {
        if (es.eigenvalues()[i]<1e-6) {
            DiagEval(i,i) = 1e-3;
        }
    }
    
    H = Evec * DiagEval * Evec.transpose();
    
}
