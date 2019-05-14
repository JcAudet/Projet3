#ifndef DOMAIN_H
#define DOMAIN_H


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
        virtual ~Domain();

        // Getters

        unsigned int GetTotNumPoint(){return totNumPoints_;}

        // Setters

        // Method

    protected:

    private:



};

}; // Namespace: diff_elec

#endif // DOMAIN_H
