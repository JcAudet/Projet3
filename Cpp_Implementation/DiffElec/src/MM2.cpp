
#include"../include/DiffElec_bits/MM2.h"
#include<iostream>

namespace DiffElec
{

MM2::MM2(Domain2D* dom, Pot* V, double dt) :
    V_(V), dom_(dom), dx_(dom->getDx(0)), dy_(dom->getDx(1)), dt_(dt)
{
    // Getting 2D potential in 1D rep
    arma::vec VVec = V_->getVectorized();

    // Computation of coefficient
    b_ = std::complex<double>(1,0) + I*HBAR*dt_*( 1/pow(dx_,2) + 1/pow(dy_,2) )/( 2*ME ) + ( I*dt_/(2*HBAR) )* VVec;
    f_ = std::complex<double>(1,0) - I*HBAR*dt_*( 1/pow(dx_,2) + 1/pow(dy_,2) )/( 2*ME ) - ( I*dt_/(2*HBAR) )* VVec;

    c_ = -I*HBAR*dt_ / (4*ME*pow(dx_,2));
    d_ = -I*HBAR*dt_ / (4*ME*pow(dy_,2));

    g_ = I*HBAR*dt_ / (4*ME*pow(dx_,2));
    k_ = I*HBAR*dt_ / (4*ME*pow(dy_,2));

    // Calculation of tot num of element in sparse
    totElem_ = 5*VVec.size() - 2 - 2*dom_->getTotNumPointsPerDim(1);

    // Creation of sparse location and values
    arma::umat pos = arma::umat(2,totElem_,arma::fill::zeros);
    arma::cx_vec val = arma::cx_vec(totElem_,arma::fill::zeros);
    arma::cx_vec val2 = arma::cx_vec(totElem_,arma::fill::zeros);

    int j=0;
    for(unsigned int i = 0; i<VVec.size(); i++)
    {
        pos(0,j) = i;
        pos(1,j) = i;
        val(j) = VVec(i);
        val2(j) = VVec(i);

        j+=1;
    }
    for(unsigned int i = dom_->getTotNumPointsPerDim(1); i<VVec.size(); i++)
    {
        pos(0,j) = i - dom_->getTotNumPointsPerDim(1);
        pos(1,j) = i;
        val(j) = c_;
        val2(j) = g_;

        j+=1;

        pos(0,j) = i;
        pos(1,j) = i - dom_->getTotNumPointsPerDim(1);
        val(j) = c_;
        val2(j) = g_;

        j+=1;
    }
    for(unsigned int i = 0; i<VVec.size() - 1; i++)
    {
        pos(0,j) = i + 1;
        pos(1,j) = i;

        if( i % dom_->getTotNumPointsPerDim(1)!=0 && i!=0 ){
            val(j) = d_;
            val2(j) = k_;
        }
        else {
            val(j) = 0;
            val2(j) = 0;
        }

        j+=1;

        pos(0,j) = i;
        pos(1,j) = i + 1;

        if( i % dom_->getTotNumPointsPerDim(1)!=0 && i!=0 ){
            val(j) = d_;
            val2(j) = k_;
        }
        else {
            val(j) = 0;
            val2(j) = 0;
        }

        j+=1;

    }

    M_ = arma::sp_cx_mat(pos,val);

    M2_ = arma::sp_cx_mat(pos,val2);

}

bool MM2::solve(arma::cx_vec* psiV)
{
    arma::cx_vec B = M2_ * (*psiV);

    arma::superlu_opts opts;
    opts.allow_ugly = true;

    bool success = spsolve(*psiV, M_, B, "superlu", opts);

    std::cout << success << std::endl;

    return success;
}

} //namespace DiffElec
