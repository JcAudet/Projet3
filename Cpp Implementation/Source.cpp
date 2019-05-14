
#include<algorithm>
#include<vector>
#include<boost/math/constants/constants.hpp>
#include<boost/numeric/ublas/vector.hpp>
#include<boost/numeric/ublas/io.hpp>


boost::numeric::ublas::vector<double> linspace(double sf,double ds) {

	boost::numeric::ublas::vector<double> S(floor(sf / ds));
	std::generate(S.begin(), S.end(), [ds]() {
		static double i = 0;
		return i += ds;
	});
	return S;
}

int main() {

	// Constantes
	double hbar = 1.054e-34;
	double m = 9.019e-31;
	double PI = boost::math::constants::pi<double>();

	// Discretisation
	double dx = 5e-12;
	double xf = 1e-9;
	double dy = dx;
	double yf = 1e-9;

	vector<double> X = linspace(xf, dx);
	vector<double> Y = linspace(yf, dy);

	// Parametres
	double sig = 8e-11;
	double lam = 5e-11;
	std::vector<double> pos = {1,2};
	double k = (2 * PI) / lam;
	std::vector<double> kv(2); kv[0] = k;


	system("pause");
	return 0;
}
