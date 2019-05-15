#ifndef MM2_H
#define MM2_H

#include<armadillo>
#include<constants.h>
#include<Potential.h>

namespace DiffElec
{

class MM2
{
    public:
        // Tructors
        MM2(Domain2D* dom, Pot* V, double dx, double dy, double dt);
        virtual ~MM2(){};

        // Getters
        arma::cx_vec getCoeff();

        // Methods



    protected:

    Domain2D* dom_;

    sp_cx_mat M_; // sp_cx_mat --> SpMat<std::complex<double>>
    sp_cx_mat M2_;

    Pot* V_;
    double dx_;
    double dy_;
    double dt_;

    complex<double> b_;
    complex<double> f_;
    complex<double> c_;
    complex<double> d_;
    complex<double> g_;
    complex<double> k_;

};

} // namespace DiffElec

#endif // MM2_H
