#include "../include/DiffElec_bits/Pot.h"


namespace DiffElec
{

Pot::Pot(Domain2D* dom) : dom_(dom)
{
    V_ = new arma::mat( dom_->getTotNumPointsPerDim(1), dom_->getTotNumPointsPerDim(0), arma::fill::zeros);
    //V_ = new arma::mat( dom_->getTotNumPointsPerDim(1), dom_->getTotNumPointsPerDim(0), arma::fill::ones);
}

} // namespace DiffElec
