#pragma once

#include"constants.h"
#include"Pot.h"
#include"Psi.h"

namespace DiffElec
{

class MM2
{
    public:
        // Tructors
        MM2(Domain2D* dom, Pot* V, double dt);
        virtual ~MM2(){};

        // Getters
        arma::cx_vec getCoeff();

        // Methods
        bool solve(arma::cx_vec* psiV);


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

