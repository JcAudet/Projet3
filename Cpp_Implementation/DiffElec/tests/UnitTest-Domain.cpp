
/*!
* /class Domain2D
* /author Jean-Christophe Gagnon-Audet
* /date 14/05/19
* /brief Unit Test for the Domain class
*/

#include<gtest/gtest.h>
// #include<diffelec>

// using namespace diff_elec;

TEST(Domain,Domain2D)
{
    Domain2D* myDom = new Domain2D(100,100,-10,10,-20,20);

    ASSERT_EQ(100,myDom->getTotNumPointsPerDim(0));
    ASSERT_EQ(100,myDom->getTotNumPointsPerDim(1));
    ASSERT_EQ(100*100,myDom->getTotNumPoints());

    ASSERT_NEAR(-10.0,myDom->getMin(0),1.0e-12);
    ASSERT_NEAR(-10.0,myDom->getMin(1),1.0e-12);
    ASSERT_NEAR(10.0,myDom->getMax(0),1.0e-12);
    ASSERT_NEAR(10.0,myDom->getMax(1),1.0e-12);

    ASSERT_DOUBLE_EQ(20.0/99.0,myDom->getDx(0));
    ASSERT_DOUBLE_EQ(40.0/99.0,myDom->getDx(1));
    ASSERT_DOUBLE_EQ(20.0*20.0/(99.0*99.0),myDom->getDiffVolumeEl());

    for*std::size_t i=0 ; i<50 ; i++)
    {
        double value_x = myDom->getMin(0) + i*myDom->getDx(0);
        double value_y = myDom->getMin(1) + i*myDom->getDx(1);

        ASSERT_EQ(value_x,myDom->getPtCoords(0,i));
        ASSERT_EQ(value_y,myDom->getPtCoords(1,i));
    }

    delte myDom;

}

GTEST_API_ int main(int argc, char **argv)
{
    printf("Running main() UnitTest-Domain.cpp. \n");
    testing::InitGoogleTest(&argc, argv);
    auto result = RUN_ALL_TESTS();

    return result;
}
