#ifndef MM2_H
#define MM2_H

#include<armadillo>

class MM2
{
    public:
        // Tructors
        MM2(arma::mat V, double dx, double dy, double dt);
        virtual ~MM2(){};

        // Getters
        arma::cx_vec getCoeff();

        // Methods



    protected:

    sp_cx_mat M; // sp_cx_mat --> SpMat<std::complex<double>>
    sp_cx_mat M2;

    complex<double> b_;
    complex<double> f_;
    complex<double> c_;
    complex<double> d_;
    complex<double> g_;
    complex<double> k_;

};

#endif // MM2_H
