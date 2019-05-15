

#include "../include/DiffElec_bits/Psi.h"

namespace DiffElec
{

Psi::Psi(Domain2D* dom, MM2* mmat, double sigmaX, double sigmaY, double kX, double kY, double posX0, double posY0, std::string nameOfFile):
    dom_(dom), nameOfFile_(nameOfFile)
{


    // Initialisation of parameters
    sigma_.push_back(sigmaX);
    sigma_.push_back(sigmaY);
    k_.push_back(kX);
    k_.push_back(kY);
    pos0_.push_back(posX0);
    pos0_.push_back(posY0);

    // Initialisation of Wave Packet
    arma::cx_rowvec tempPsiX = exp(-arma::pow( dom_->getCoords(0)-pos0_[0] ,2 ) / pow( sigma_[0],2 ) ) * exp(I * k_[0] * dom_->getCoords(0) );
    arma::cx_colvec tempPsiY = exp(-arma::pow( dom_->getCoords(1)-pos0_[1] ,2 ) / pow( sigma_[1],2 ) ) * exp(I * k_[1] * dom_->getCoords(1) );

    arma::cx_mat psiMat = tempPsiY * tempPsiX;
    psiV_ = arma::vectorise(psiMat);

}

void Psi::iterate()
{

    mmat_->solve(&psiV_); // need output
    save();

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
