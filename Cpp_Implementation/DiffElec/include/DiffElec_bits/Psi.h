#pragma once

#include<string>
#include"Domain.h"

namespace DiffElec
{

class Psi
{
    public:

        // tructors
        Psi(Domain2D* dom, double sigmaX, double sigmaY, double kX, double kY, double posX0, double posY0);
        virtual ~Psi(){}

        // Getters
        std::vector<double> getSigma(){return sigma_;};
        std::vector<double> getK(){return k_;};
        std::vector<double> getpos(){return pos0_;};

        // Method
        void iterate(const std::string nameOfFile);
        void save(const std::string nameOfFile);
        double normalizeVec[...]

    protected:

        Domain2D* dom_;

        arma::vec psiI_;
        arma::vec psiF_;

        std::vector<double> sigma_;
        std::vector<double> k_;
        std::vector<double> pos0_;

};

}; // Namespace DiffElec

