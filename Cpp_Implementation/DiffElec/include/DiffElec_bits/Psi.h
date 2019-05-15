#pragma once

#include"MM2.h"

namespace DiffElec
{

class Psi
{
    public:

        // tructors
        Psi(Domain2D* dom, MM2* mmat, double sigmaX, double sigmaY, double kX, double kY, double posX0, double posY0);
        virtual ~Psi(){}

        // Getters
        std::vector<double> getSigma(){return sigma_;};
        std::vector<double> getK(){return k_;};
        std::vector<double> getpos(){return pos0_;};

        // Method
        void iterate(const std::string nameOfFile);
        void save(const std::string nameOfFile);
        double normalizeVec();
        double trapeze2D();

    protected:

        Domain2D* dom_;

        MM2* mmat_;

        arma::vec psiV_;
        arma::mat psiM_;

        std::vector<double> sigma_;
        std::vector<double> k_;
        std::vector<double> pos0_;

};

}; // Namespace DiffElec

