#include "GaussIntegration.h"

void GaussIntegration::init(int n) {
	w.resize(n);
	var.resize(n);

	switch (n)
	{
	case 2:
		w[0] = w[1] = 1.0;
		var[0] = -0.577350269189625;
		var[1] = -var[0];
		break;
	case 3:
		w[0] = 0.555555556;
		w[1] = 0.888888889;
		w[2] = w[1];
		var[0] = -0.774596669;
		var[1] = 0.0;
		var[2] = -var[0];
		break;
	case 5:
		w[0] = 0.2369268851;
		w[1] = 0.4786286705;
		w[2] = 0.5688888888;
		w[3] = w[1];
		w[4] = w[0];
		var[0] = -0.9061798459;
		var[1] = -0.5384693101;
		var[2] = 0.0;
		var[3] = -var[1];
		var[4] = -var[0];
		break;
	default:
		throw "Gauss Integration not implemented for this order N.";
	}
}


double GaussIntegration::nPointsGauss(double a, double b, std::function<double(double)> f) {
	double result = 0.0;
	for (int i = 0; i < w.size(); i++) {
		result += w[i] * f(varChange(var[i], a, b));
	}
	return result *= (b - a) / 2.0;
}

double GaussIntegration::nPointsGauss(double a, double b, double c, double d, std::function<double(double, double)> f) {
	double res = 0;
	for (int i = 0; i < w.size(); i++) {
		for (int j = 0; j < w.size(); j++) {
			res += w[i] * w[j] * f(varChange(var[i], a, b), varChange(var[j], c, d));
		}
	}
	res *= ((b - a) / 2.0) * ((d - c) / 2.0);
	return res;
}