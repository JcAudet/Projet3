

#include "../include/DiffElec_bits/Psi.h"

#include<iostream>

namespace DiffElec
{

Psi::Psi(Domain2D* dom, MM2* mmat, double sigmaX, double sigmaY, double kX, double kY, double posX0, double posY0, std::string nameOfFile):
    dom_(dom), mmat_(mmat),nameOfFile_(nameOfFile)
{

    // Initialisation of parameters
    sigma_.push_back(sigmaX);
    sigma_.push_back(sigmaY);
    k_.push_back(kX);
    k_.push_back(kY);
    pos0_.push_back(posX0);
    pos0_.push_back(posY0);

    // Initialisation of Wave Packet
    //auto tempPsiX = exp(- arma::pow( dom_->getCoords(0)-pos0_[0] , 2 ) / pow( sigma_[0],2 ) ) % exp(I * k_[0] * dom_->getCoords(0) );
    //auto tempPsiY = exp(- arma::pow( dom_->getCoords(1)-pos0_[1] , 2 ) / pow( sigma_[1],2 ) ) % exp(I * k_[1] * dom_->getCoords(1) );

    arma::cx_rowvec tPX = arma::cx_rowvec(dom_->getTotNumPointsPerDim(0),arma::fill::zeros);
    for(unsigned int i = 0; i<dom_->getTotNumPointsPerDim(0);i++)
    {
        tPX(i) = exp(-pow( dom_->getPtCoords(0,i)-pos0_[0] , 2 ) / pow( sigma_[0],2 ) ) * exp(I * k_[0] * dom_->getPtCoords(0,i) );
    }
    arma::cx_vec tPY = arma::cx_vec(dom_->getTotNumPointsPerDim(1),arma::fill::zeros);
    for(unsigned int i = 0; i<dom_->getTotNumPointsPerDim(1);i++)
    {
        tPY(i) = exp(-pow( dom_->getPtCoords(1,i)-pos0_[1] , 2 ) / pow(sigma_[1],2) ) * exp(I * k_[1] * dom_->getPtCoords(1,i) );
    }

    arma::cx_mat psiMat = tPY * tPX;

    // Normalisation and storing
    norm_ = normaliseMat(psiMat);

    psiV_ = arma::vectorise(psiMat);

}

void Psi::iterate()
{

    std::cout << arma::size(psiV_) << std::endl;
    bool success = mmat_->solve(&psiV_);
    std::cout << arma::size(psiV_) << std::endl;

    save();

}

double Psi::normaliseMat(arma::cx_mat& psi)
{
    double norm = dom_->getDiffVolumeEl() * arma::sum(arma::sum(arma::pow(arma::abs(psi),2)));
    psi /= norm;
    return dom_->getDiffVolumeEl() * arma::sum(arma::sum(arma::pow(arma::abs(psi),2)));
}

void Psi::save()
{
    std::string nameOfFileDat = nameOfFile_ + ".dat";
    std::ofstream myFile;
    myFile.open(nameOfFileDat);
    myFile << std::scientific;
    myFile.precision(15);
    psiV_.raw_print(myFile);
    myFile.close();
};


}; // namespace DiffElec
