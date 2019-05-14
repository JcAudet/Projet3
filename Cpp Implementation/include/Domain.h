#ifndef DOMAIN_H
#define DOMAIN_H

#include<vector>
#include<string>
#include<armadillo>

namespace diff_elec
{

/*!
* /class Domain
* /author Jean-Christophe Gagnon-Audet
* /date 14/05/19
* /brief This class defines the finite domain for calculation
*/

class Domain
{
    public:

        // Virtual Destructor
        virtual ~Domain(){};

        // Getters
        unsigned int getTotNumPoint(){return totNumPoints_;}
        unsigned int getTotNumPointsPerDim(unsigned int dim){return numPointsPerDim_[dim];}
        double getMin(unsigned int dim){return minX_[dim];}
        double getMax(unsigned int dim){return maxX_[dim];}
        double getDx(unsgined int dim){return dx_[dim];}
        double getDiffVolumeEl(){return diffVolumeEl_;}
        double getPtCoords(unsigned int dim, unsigned int pt){return ptCoords_[dim][pt];}
        unsigned int getDim(){return dim_;}

        // Method
        void save(const std::string nameOfFile)

    protected:

        // Default Constructor
        Domain(int numThreads, unsigned int dim)

        std::vector<arma::vec> ptCoords_;
        std::vector<unsgined int> numPointsPerDim_;
        unsigned int totNumPoints_;
        unsigned int dim_;
        std::vector<double> minX_;
        std::vector<double> maxX_;
        std::vector<double> dx_;
        double diffVolumeEl_;

};

/*!
* /class Domain2D
* /author Jean-Christophe Gagnon-Audet
* /date 14/05/19
* /brief This class defines the finite 2D domain for calculation
*/

class Domain2D : public Domain
{
    public:

        //Default constructor
        Domain2D(unsigned int numPointX1,unsigned int numPointX2,
                double x1Min, double x1Max,
                double x2Min, double x2Min,
                int numThreads);

        // Default Destructor
        virtual ~Domain2D(){};
}


}; // Namespace: diff_elec

#endif // DOMAIN_H