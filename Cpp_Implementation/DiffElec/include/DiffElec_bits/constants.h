
#pragma once

#include<complex>
#include<boost/units/systems/si/codata/universal_constants.hpp>
#include<boost/units/systems/si/codata/electron_constants.hpp>
#include <boost/math/constants/constants.hpp>

namespace DiffElec
{

const std::complex<double> I = std::complex<double>(0,1);
//const double HBAR = boost::units::si::constants::codata::hbar.value();
//const double ME = boost::units::si::constants::codata::m_e.value();
const double HBAR = 1.0545718e-34;
const double ME = 9.10938e-31;
const double PI = boost::math::constants::pi<double>();

} // namespace DiffElec
