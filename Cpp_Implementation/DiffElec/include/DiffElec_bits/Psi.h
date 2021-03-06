#pragma once

#include"Domain.h"
#include"constants.h"
#include"MM2.h"

namespace DiffElec
{
class MM2;
class Psi
{
    public:

        // tructors
        Psi(Domain2D* dom, MM2* mmat, double sigmaX, double sigmaY, double kX, double kY, double posX0, double posY0, std::string nameOfFile);
        virtual ~Psi(){}

        // Getters
        std::vector<double> getSigma(){return sigma_;};
        std::vector<double> getK(){return k_;};
        std::vector<double> getpos(){return pos0_;};
        double getNorm(){return norm_;}

        // Method
        void iterate();
        void save();
        double trapezoid(arma::cx_mat psi);
        double normaliseMat(arma::cx_mat& psi);
        arma::cx_mat vec2mat(arma::cx_vec psi);

    protected:

        Domain2D* dom_;
        MM2* mmat_;

        std::string nameOfFile_;

        arma::cx_vec psiV_;

        double norm_;

        std::vector<double> sigma_;
        std::vector<double> k_;
        std::vector<double> pos0_;

};

}; // Namespace DiffElec

