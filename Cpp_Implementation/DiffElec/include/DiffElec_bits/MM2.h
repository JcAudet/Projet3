#pragma once

#include<complex>
#include"constants.h"
#include"Pot.h"

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

    Pot* V_;
    Domain2D* dom_;
    double dx_;
    double dy_;
    double dt_;

    arma::cx_vec b_;
    arma::cx_vec f_;
    std::complex<double> c_;
    std::complex<double> d_;
    std::complex<double> g_;
    std::complex<double> k_;

    unsigned int totElem_;

    arma::sp_cx_mat M_; // sp_cx_mat --> SpMat<std::complex<double>>
    arma::sp_cx_mat M2_;

};

}; // namespace DiffElec

