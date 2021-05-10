#include <string>
#include <vector>

void vectorSubtract(std::vector<double>& a, std::vector<double>& b, std::vector<double>& result) {
	for (int i = 0; i < result.size(); i++) result[i] = a[i] - b[i];
}

void vectorAdd(std::vector<double>& a, std::vector<double>& b, std::vector<double>& result) {
	for (int i = 0; i < result.size(); i++) result[i] = a[i] + b[i];
}

double scalarProduct(std::vector<double>& a, std::vector<double>& b) {
	double sum = 0;
	for (int i = 0; i < a.size(); i++)	sum += a[i] * b[i];
	return sum;
}

void vectorMultOnConst(std::vector<double>& a, double constant, std::vector<double>& result) {
	for (int i = 0; i < result.size(); i++) result[i] = a[i] * constant;
}

void vectorAssigment(std::vector<double>& a, std::vector<double>& b) {
	for (int i = 0; i < a.size(); i++) a[i] = b[i];
}

double calcNorm(std::vector<double>& vector) {
	double sum = 0;
	for (int i = 0; i < vector.size(); i++) sum += pow(vector[i], 2);
	return sqrt(sum);
}