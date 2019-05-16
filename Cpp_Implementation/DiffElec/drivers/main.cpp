#include "../include/diffelec"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

int main()
{

    // Parsing parameters
    boost::property_tree::ptree pt;
    boost::property_tree::ini_parser::read_ini("parameters.ini", pt);

    // Constructing
    DiffElec::Domain2D* dom = new DiffElec::Domain2D(pt.get<unsigned int>("Domain.nbX"), pt.get<unsigned int>("Domain.nbY"),
                    pt.get<double>("Domain.xI"), pt.get<double>("Domain.xF"),
                    pt.get<double>("Domain.yI"), pt.get<double>("Domain.yF"));

    DiffElec::Pot* V = new DiffElec::Pot(dom);
    DiffElec::MM2* mmat = new DiffElec::MM2(dom, V, pt.get<double>("Time.dt"));

    DiffElec::Psi* psi = new DiffElec::Psi(dom, mmat, pt.get<double>("Wave_Packet.sigmaX"), pt.get<double>("Wave_Packet.sigmaY"),
                        pt.get<double>("Wave_Packet.kX"), pt.get<double>("Wave_Packet.kY"),
                        pt.get<double>("Wave_Packet.posX0"), pt.get<double>("Wave_Packet.posX0"),
                        pt.get<std::string>("Other.nameOfFile"));

    // Evolve
    for(unsigned int i=0; i< pt.get<unsigned int>("Time.timeSteps"); i++)
    {
        psi->iterate();
    }

    delete dom;
    delete V;
    delete mmat;
    delete psi;






}
