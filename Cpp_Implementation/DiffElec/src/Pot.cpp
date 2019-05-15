#include "../include/DiffElec_bits/Pot.h"


namespace DiffElec
{

Pot::Pot(Domain2D* dom) : dom_(dom)
{
    arma::mat V_( dom_->getTotNumPointsPerDim(1), dom_->getTotNumPointsPerDim(0), arma::fill::zeros);
}

} // namespace DiffElec
