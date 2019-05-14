#include "Domain.h"

namespace diff_elec
{

Domain::Domain(int numThreads, unsigned int dim) : dim_(dim)
{
    ptCoords_.resize(dim_);
    dx.resize(dim_);
}

void Domain::save(std::string nameOfFile)
{
    for(std::size_t i=0; i<dim_; i++)
    {
        std::string nameOfFileDim = nameOfFile + "_dim" + std::to_string(i) + ".dat";
        std::ofstream myFile;
        myFile.open(nameOfFileDim);
        myFile << std::scientific;
        myFile.precision(15);
        ptCoords_[i].raw_print(myFile);
        myFile.close();
    }
}

Domain2D::Domain2D(unsigned int numPointX1,unsigned int numPointX2,
                double x1Min, double x1Max,
                double x2Min, double x2Min,
                int numThreads=1) : Domain(numThreads,2)
{
    numPointsPerDim_.push_back(numPointX1);
    numPointsPerDim_.push_back(numPointX2);
    minX_.push_back(x1Min);
    minX_.push_back(x2Min);
    maxX_.push_back(x1Max);
    maxX_.push_back(x2Max);
    totNumPoints_=1;
    diffVolumeEl_=1.0;
    for(std::size_t i=0 ; i<dim_ ; i++)
    {
        // Compute total num of points
        totNumPoints_*= numPointsPerDim_[i];

        // Compute distance between each point
        if(numPointsPerDim_[i] == 1) dx[i]=1.0;
        else dx[i] = std::abs(maxX_[i]-minX_[i]) / double(numPointsPerDim_[i]-1);

        // Compute the diffential volume elements
        diffVolumeEl_ *= dx_[i];

        // Compute the position of points
        ptCoords_[i].zeros(numPointsPerDim_[i]);
        for(std::size_t nx=0;i<numPointsPerDim_[i]; nx++)
        {
            ptCoords_[i][nx] = minX_[i] + nx*dx_[i];
        }
    }

}

} // namespace diff_elec
