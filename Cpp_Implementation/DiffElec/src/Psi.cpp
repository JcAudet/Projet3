

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

    std::cout << norm_ << std::endl;

    psiV_ = arma::vectorise(psiMat);

}

void Psi::iterate()
{

    bool success = mmat_->solve(&psiV_);

    //save();

}

double Psi::normaliseMat(arma::cx_mat& psi)
{
    double norm = dom_->getDiffVolumeEl() * arma::sum(arma::sum(arma::pow(arma::abs(psi),2)));
    psi = psi / norm;
    return dom_->getDiffVolumeEl() * arma::sum(arma::sum(arma::pow(arma::abs(psi),2)));
}

void Psi::save()
{

    arma::cx_mat fP = vec2mat(psiV_);

    std::string nameOfFileDat = nameOfFile_ + ".dat";
    std::ofstream myFile;
    myFile.open(nameOfFileDat,std::ios_base::app);
    myFile << std::scientific;
    myFile.precision(15);
    arma::pow( arma::abs( fP ) , 2 ).raw_print(myFile);
    //myFile << "////// \n";
    myFile.close();
};

arma::cx_mat Psi::vec2mat(arma::cx_vec psi)
{
    std::cout << dom_->getTotNumPointsPerDim(1) << std::endl;

    size_t rows = std::ceil(psi.n_elem / double(dom_->getTotNumPointsPerDim(1)));

    std::cout << rows << std::endl;

    arma::cx_mat m = reshape(psi, dom_->getTotNumPointsPerDim(1), rows);

    return m;

}

}; // namespace DiffElec
