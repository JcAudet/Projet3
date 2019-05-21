#include "../include/DiffElec_bits/Pot.h"


namespace DiffElec
{

Pot::Pot(Domain2D* dom) : dom_(dom)
{
    V_ = new arma::mat( dom_->getTotNumPointsPerDim(1), dom_->getTotNumPointsPerDim(0), arma::fill::zeros);
}


Pot1Ft::Pot1Ft(Domain2D* dom, double V, double xm, double h, double ym, double a): Pot(dom)
{

    V_ = new arma::mat( dom_->getTotNumPointsPerDim(1), dom_->getTotNumPointsPerDim(0), arma::fill::zeros);

    for(unsigned int i=0; i < dom_->getTotNumPointsPerDim(0); i++){
        if(dom_->getPtCoords(0,i) > xm - h/2 && dom_->getPtCoords(0,i) < xm + h/2){
            for(unsigned int j=0; j < dom_->getTotNumPointsPerDim(1); j++){
                if(dom_->getPtCoords(1,j) > ym-a/2 && dom_->getPtCoords(1,j) < ym+a/2){
                    (*V_)(j,i) = V;
                }
            }
        }
    }

}

} // namespace DiffElec
