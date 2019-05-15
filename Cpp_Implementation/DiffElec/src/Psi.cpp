

#include "../include/DiffElec_bits/Psi.h"

namespace DiffElec
{

Psi::Psi(Domain2D* dom, double sigmaX, double sigmaY, double kX, double kY, double posX0, double posY0):
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
    //arma::rowvec tempPsiX = exp(-pow( dom_->getCoords(0)-pos0_[0] ,2 ) / pow( sigma_[0],2 ) ) * exp(std::complex<double>(0,1) * k_[0] * dom_->getCoords(0));
    //arma::colvec tempPsiY = exp(-pow( dom_->getCoords(1)-pos0_[1] ,2 ) / pow( sigma_[1],2 ) ) * exp(std::complex<double>(0,1) * k_[1] * dom_->getCoords(1));

    // psiIni = tempPsiY * tempPsiX (OUTER PRODUCT)
    // psiI = arma::vectorize(psiIni)


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
