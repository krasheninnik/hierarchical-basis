#pragma once
#include <vector>
#include <functional>

class GaussIntegration {
public:
	void init(int n);
	//x from [a,b], n - Gauss points number
	double nPointsGauss(double a, double b, std::function<double(double)>);
	//x from [a,b], y form [c,d], n - Gauss points number
	double nPointsGauss(double a, double b, double c, double d, std::function<double(double, double)>);
	double varChange(double var, double a, double b);

private:
	//double varChange(double var, double a, double b);

	std::vector<double> w;
	std::vector<double> var;
};