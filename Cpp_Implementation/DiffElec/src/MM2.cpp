
#include<../include/DiffElec_bits/MM2.h>

namespace DiffElec
{

MM2::MM2(Domain2D* dom, Pot* V, double dx, double dy, double dt) :
    V_(V), dom_(dom), dx_(dx), dy_(dy), dt_(dt)
{
    // Getting 2D potential in 1D rep
    arma::vec VVec = V_->getVectorized();

    // Computation of coefficient
    b_ = 1 + I*HBAR*dt_*( 1/dx_^2 + 1/dy_^2 )/( 2*ME ) + ( I*dt_/(2*HBAR) )* VVec;
    f_ = 1 - I*HBAR*dt_*( 1/dx_^2 + 1/dy_^2 )/( 2*ME ) - ( I*dt_/(2*HBAR) )* VVec;

    c_ = -I*HBAR*dt_ / (4*ME*dx_^2);
    d_ = -I*HBAR*dt_ / (4*ME*dy_^2);

    g_ = I*HBAR*dt_ / (4*ME*dx_^2);
    k_ = I*HBAR*dt_ / (4*ME*dy_^2);

    // Calculation of tot num of element in sparse
    unsigned int totElem = 5*VVec.size() - 2 - 2*dom_->getTotNumPointsPerDim(1);

    // Creation of sparse location and values
    //umat pos = umat(2,totElem);
    //cx_vec val = cx_vec(totElem);
    umat pos;
    cx_vec val,val2;
    for(unsigned int i = 0; i<VVec.size(); i++)
    {
        pos << i << i << endr;
        val << VVec[i] << endr;
        val2 << VVec[i] << endr;
    }
    for(unsigned int i = dom_->getTotNumPointsPerDim(1); i<VVec.size(); i++)
    {
        pos << i - dom_->getTotNumPointsPerDim(1) << i << endr;
        pos << i << i - dom_->getTotNumPointsPerDim(1) << endr;
        val << c_ << endr << c << endr;
        val2 << g_ << endr << c << endr;
    }
    for(unsigned int i = 0; i<VVec.size() - 1; i++)
    {
        pos << i + 1 << i << endr;
        pos << i << i + 1 << endr;
        if( i % dom_->getTotNumPointsPerDim(1)!=0 && i!=0 ){
            val << d_ << endr << d_<< endr;
            val2 << k_ << endr << k_ << endr;
        }
        else {
            val << 0 << endr << 0 << endr;
            val2 << 0 << endr << 0 << endr;
            }
    }

    M_ = sp_cx_mat(pos,val);
    M2_ = sp_cx_mat(pos,val2);

}

MM2::~MM2()
{
    //dtor
}

} //namespace DiffElec
