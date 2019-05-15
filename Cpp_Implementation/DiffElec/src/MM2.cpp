
#include"../include/DiffElec_bits/MM2.h"

namespace DiffElec
{

MM2::MM2(Domain2D* dom, Pot* V, double dx, double dy, double dt) :
    V_(V), dom_(dom), dx_(dx), dy_(dy), dt_(dt)
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
    arma::umat pos;
    arma::cx_vec val,val2;
    for(unsigned int i = 0; i<VVec.size(); i++)
    {
        pos << i << i << arma::endr;
        val << VVec[i] << arma::endr;
        val2 << VVec[i] << arma::endr;
    }
    for(unsigned int i = dom_->getTotNumPointsPerDim(1); i<VVec.size(); i++)
    {
        pos << i - dom_->getTotNumPointsPerDim(1) << i << arma::endr;
        pos << i << i - dom_->getTotNumPointsPerDim(1) << arma::endr;
        val << c_ << arma::endr << c_ << arma::endr;
        val2 << g_ << arma::endr << g_ << arma::endr;
    }
    for(unsigned int i = 0; i<VVec.size() - 1; i++)
    {
        pos << i + 1 << i << arma::endr;
        pos << i << i + 1 << arma::endr;
        if( i % dom_->getTotNumPointsPerDim(1)!=0 && i!=0 ){
            val << d_ << arma::endr << d_<< arma::endr;
            val2 << k_ << arma::endr << k_ << arma::endr;
        }
        else {
            val << 0 << arma::endr << 0 << arma::endr;
            val2 << 0 << arma::endr << 0 << arma::endr;
            }
    }

    M_ = arma::sp_cx_mat(pos,val);
    M2_ = arma::sp_cx_mat(pos,val2);

}

void MM2::solve(arma::cx_vec* psiV)
{
    //LINEAR SOLVE!
}

} //namespace DiffElec
