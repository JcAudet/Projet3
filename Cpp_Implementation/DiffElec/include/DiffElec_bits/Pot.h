#pragma once
#pragma message("Pot Included")

#include"Domain.h"

namespace DiffElec
{

class Pot
{
    public:

        // Tructors
        Pot(Domain2D* dom);
        virtual ~Pot(){}

        // Getters
        arma::mat getV(){return V_;}
        arma::vec getVectorized(){
            return arma::vectorise(V_,0);
        }

    protected:

        Domain2D* dom_;
        arma::mat V_;

};

}; // namespace DiffElec
