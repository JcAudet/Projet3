
/*!
* /class Domain2D
* /author Jean-Christophe Gagnon-Audet
* /date 14/05/19
* /brief Unit Test for the Pot class
*/

#include<gtest/gtest.h>
#include"../include/diffelec"

#include<iostream>

using namespace DiffElec;

TEST(Pot,Pot)
{
    Domain2D* mydom = new Domain2D(100,100,-10.0,10.0,-20.0,20.0,1);

    Pot* mypot = new Pot(mydom);

    ASSERT_EQ(100*100, mypot -> getVectorized().size());

    delete mydom;
    delete mypot;

}

GTEST_API_ int main(int argc, char **argv)
{
    printf("Running main() UnitTest-Pot.cpp. \n");
    testing::InitGoogleTest(&argc, argv);
    auto result = RUN_ALL_TESTS();

    return result;
}
