

#include "../include/DiffElec_bits/Psi.h"

namespace DiffElec
{

Psi::Psi(Domain2D* dom, MM2* mmat, double sigmaX, double sigmaY, double kX, double kY, double posX0, double posY0):
    dom_(dom)
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

    psiM_ = tempPsiY * tempPsiX;
    psiV_ = arma::vectorise(psiM_);


}

void Psi::iterate(const std::string nameOfFile)
{
    psiI_.clear();

    this->save(nameOfFile);

    psiI_=psiF_;


}

void Psi::save(std::string nameOfFile)
{
    std::string nameOfFileDat = nameOfFile + ".dat";
    std::ofstream myFile;
    myFile.open(nameOfFileDat);
    myFile << std::scientific;
    myFile.precision(15);
    psiF_.raw_print(myFile);
    myFile.close();
};


}; // namespace DiffElec
