

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

    // Creation of initiale state
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
    // Outer product
    arma::cx_mat psiMat = tPY * tPX;

    // Normalisation and storing
    norm_ = normaliseMat(psiMat);
    psiV_ = arma::vectorise(psiMat);

}

void Psi::iterate()
{

    arma::cx_vec temp = mmat_->solve(psiV_);

    psiV_ = temp;

    //save();

}

double Psi::normaliseMat(arma::cx_mat& psi)
{
    double norm = trapezoid(psi);
    psi = psi / std::sqrt(norm);
    return trapezoid(psi);
}

double Psi::trapezoid(arma::cx_mat Psi)
{

    arma::mat psi = arma::pow( arma::abs(Psi) , 2);

    double sumABCD = psi(0,0) + psi(0, dom_->getTotNumPointsPerDim(0)-1) + psi(dom_->getTotNumPointsPerDim(1)-1,0) + psi(dom_->getTotNumPointsPerDim(1)-1,dom_->getTotNumPointsPerDim(0)-1);
    double sumA = arma::sum ( psi.row(0) ) - psi(0,0) - psi(0,dom_->getTotNumPointsPerDim(0)-1);
    double sumB = arma::sum ( psi.row(dom_->getTotNumPointsPerDim(1)-1) ) - psi(dom_->getTotNumPointsPerDim(1)-1,0) - psi(dom_->getTotNumPointsPerDim(1)-1,dom_->getTotNumPointsPerDim(0)-1);
    double sumC = arma::sum ( psi.col(0) ) - psi(dom_->getTotNumPointsPerDim(1)-1,0) - psi(0,0);
    double sumD = arma::sum ( psi.col(dom_->getTotNumPointsPerDim(0)-1) ) - psi(0,dom_->getTotNumPointsPerDim(0)-1) - psi(dom_->getTotNumPointsPerDim(1)-1,dom_->getTotNumPointsPerDim(0)-1);
    double sumF = arma::sum( arma::sum ( psi.submat(1, 1, dom_->getTotNumPointsPerDim(1)-2, dom_->getTotNumPointsPerDim(0)-2) ) );

    return ( dom_->getDiffVolumeEl() / 4 ) * ( sumABCD + 2*(sumA + sumB + sumC + sumD) + 4*sumF );

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

    size_t rows = std::ceil(psi.n_elem / double(dom_->getTotNumPointsPerDim(1)));
    arma::cx_mat m = reshape(psi, dom_->getTotNumPointsPerDim(1), rows);

    return m;

}

}; // namespace DiffElec
