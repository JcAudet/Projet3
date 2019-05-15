#include "../include/DiffElec_bits/Potential.h"

Potential::Potential(Domain2D* dom) : dom_(dom)
{
    V_ = mat(dom_->getTotNumPerDim(0),dom_->getTotNumPerDim(1),fill::fill_zeros);
}

Potential::~Potential()
{
    //dtor
}
