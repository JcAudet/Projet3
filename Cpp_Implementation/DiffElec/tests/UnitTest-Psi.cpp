
/*!
* /class Domain2D
* /author Jean-Christophe Gagnon-Audet
* /date 14/05/19
* /brief Unit Test for the Psi class
*/

#include<gtest/gtest.h>
#include"../include/diffelec"

using namespace DiffElec;

TEST(Psi,Psi)
{
    Domain2D* mydom = new Domain2D(100,100,-10.0,10.0,-20.0,20.0,1);
    Pot* mypot = new Pot(mydom);
    MM2* mymm2 = new MM2(mydom,mypot,1);

    Psi* mypsi = new Psi(mydom, mymm2, 1, 1, 1, 1, 1, 1, "test");

    delete mydom;
    delete mypot;
    delete mymm2;
    delete mypsi;

}

GTEST_API_ int main(int argc, char **argv)
{
    printf("Running main() UnitTest-Psi.cpp. \n");
    testing::InitGoogleTest(&argc, argv);
    auto result = RUN_ALL_TESTS();

    return result;
}
