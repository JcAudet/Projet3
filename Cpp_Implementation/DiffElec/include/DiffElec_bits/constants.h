
#pragma once

#include<boost/units/systems/si/codata/universal_constants.hpp>
#include<boost/units/systems/si/codata/electron_constants.hpp>
#include <boost/math/constants/constants.hpp>

namespace DiffElec
{

using namespace boost::units::si::constants::codata;

const complex<double> I = complex<double>(0,1);
const double HBAR = value(hbar);
const double ME = value(m_e);
const double PI = boost::math::constants::pi<double>();

} // namespace DiffElec
