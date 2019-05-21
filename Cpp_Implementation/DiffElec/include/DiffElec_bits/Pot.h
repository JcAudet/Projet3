#pragma once

#include"Domain.h"

namespace DiffElec
{

class Pot
{
    public:

        // Tructors
        Pot(Domain2D* dom);
        virtual ~Pot(){
            delete V_;
        }

        // Getters
        arma::mat* getV(){return V_;}
        arma::vec getVectorized(){
            return arma::vectorise(*V_,0);
        }

    protected:

        Domain2D* dom_;
        arma::mat* V_;

};

class Pot1Ft : public Pot
{
    public:

        //Tructors
        Pot1Ft(Domain2D* dom, double V, double xm, double h, double ym, double a);
        ~Pot1Ft(){
            delete V_;
        };

};

}; // namespace DiffElec
