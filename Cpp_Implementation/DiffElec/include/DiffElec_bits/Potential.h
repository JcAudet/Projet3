#ifndef POTENTIAL_H
#define POTENTIAL_H


class Pot
{
    public:

        // Tructors
        Pot(Domain2D dom);
        virtual ~Pot();

        // Getters
        arma::mat getV(){return V;}
        arma::vec getVectorized(){
            return vectorise(V,0);
        };

    protected:

        Domain2D* dom_;
        arma::mat V_;

};


#endif // POTENTIAL_H
