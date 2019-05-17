
/*!
* /class Domain2D
* /author Jean-Christophe Gagnon-Audet
* /date 14/05/19
* /brief Unit Test for the MM2 class
*/

#include<gtest/gtest.h>
#include"../include/diffelec"

using namespace DiffElec;

TEST(MM2,MM2)
{
    Domain2D* mydom = new Domain2D(100,100,-10.0,10.0,-20.0,20.0,1);
    Pot* mypot = new Pot(mydom);
    MM2* mymm2 = new MM2(mydom,mypot,1.0);

    delete mydom;
    delete mypot;
    delete mymm2;

}

GTEST_API_ int main(int argc, char **argv)
{
    printf("Running main() UnitTest-MM2.cpp. \n");
    testing::InitGoogleTest(&argc, argv);
    auto result = RUN_ALL_TESTS();

    return result;
}
