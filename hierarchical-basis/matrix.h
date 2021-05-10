#pragma once
#include <vector>

class SparseMatrix {
	using LocalMatrix = std::vector<std::vector<double>>;
public:
	SparseMatrix();
	void init(const int elemNum, const int max_iter, const double eps);
	void init(const int dim, std::vector<int>& _ig, std::vector<int>& _jg);

	void addLocalMatrix(int elemNum, LocalMatrix M);
	int  getDim();
	int  getAmountElems();

	void addElem(int row, int column, double value);
	void addDiagElem(int n, double value);
	void fillGGU();
	void setFirstBoundaryCondition(int n);
	void resetMatrix();

	void getIG(int* source);
	void getJG(int* source);
	void getDI(double* source);
	void getGGU(double* source);
	void getGGL(double* source);
private:
	std::vector<int> ig;		// индексы элементов - начало строки
	std::vector<int> jg;		// индексы столбцов эленментов
	std::vector<double> di;		// ƒиагональные элементы
	std::vector<double> ggu;	// Ёлементы верхнего треугольника
	std::vector<double> ggl;	// Ёлементы нижнего треугольника

	int n;						// –азмерность
};