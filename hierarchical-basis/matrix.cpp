#include "matrix.h"

void SparseMatrix::getIG(int* source) {	for (int i = 0; i < ig.size(); i++) source[i] = ig[i];}

void SparseMatrix::getJG(int* source) {	for (int i = 0; i < jg.size(); i++) source[i] = jg[i];}

void SparseMatrix::getDI(double* source) {	for (int i = 0; i < di.size(); i++) source[i] = di[i];}

void SparseMatrix::getGGU(double* source) { for (int i = 0; i < ggu.size(); i++) source[i] = ggu[i]; }

void SparseMatrix::getGGL(double* source) { for (int i = 0; i < ggl.size(); i++) source[i] = ggl[i]; }

#pragma region vector_operations
inline void vectorAdd(std::vector<double>& a, std::vector<double>& b, std::vector<double>& result);
inline void vectorSubtract(std::vector<double>& a, std::vector<double>& b, std::vector<double>& result);
inline void vectorAssigment(std::vector<double>& a, std::vector<double>& b);
inline void vectorMultOnConst(std::vector<double>& a, double constant, std::vector<double>& result);
inline double scalarProduct(std::vector<double>& a, std::vector<double>& b);
inline double calcNorm(std::vector<double>& vector);
#pragma endregion

SparseMatrix::SparseMatrix() = default;
int SparseMatrix::getDim() { return n; }
int SparseMatrix::getAmountElems() { return ggl.size(); }

void SparseMatrix::resetMatrix() {
	for (auto &el : ggl) el = 0;
	for (auto &el : ggu) el = 0;
	for (auto &el : di) el = 0;
}

void SparseMatrix::init(const int dim, std::vector<int>& _ig, std::vector<int>& _jg) {
	n = dim;
	di = std::vector<double>(n);
	ig = _ig;
	jg = _jg;
	ggu = std::vector<double>(jg.size());
	ggl = std::vector<double>(jg.size());
}

void SparseMatrix::addDiagElem(int n, double value) {
	di[n] += value;
}

void SparseMatrix::addElem(int row, int column, double value) {
	int index = ig[row];
	int maxIndex = ig[row + 1];

	for (; index < maxIndex; index++) {
		if (jg[index] == column) {
			ggl[index] += value;
			return;
		}
	}

	throw "Not found column";
}

void SparseMatrix::fillGGU() {
	ggu = ggl;
}

void SparseMatrix::addLocalMatrix(int elemNum, LocalMatrix M) {
	const int diagPos = elemNum * 2;
	di[diagPos] += M[0][0];
	di[diagPos + 1] += M[1][1];
	di[diagPos + 2] += M[2][2];
	di[diagPos + 3] += M[3][3];

	const int elemPos = elemNum * 5;
	ggl[elemPos] += M[1][0];
	ggu[elemPos] += M[0][1];

	ggl[elemPos + 1] += M[2][0];
	ggu[elemPos + 1] += M[0][2];

	ggl[elemPos + 2] += M[2][1];
	ggu[elemPos + 2] += M[1][2];

	ggl[elemPos + 3] += M[3][0];
	ggu[elemPos + 3] += M[0][3];

	ggl[elemPos + 4] += M[3][1];
	ggu[elemPos + 4] += M[1][3];

	ggl[elemPos + 5] += M[3][2];
	ggu[elemPos + 5] += M[2][3];
}

void SparseMatrix::setFirstBoundaryCondition(int ind) {
	// Установка диагонального элемента = 1
	di[ind] = 1;

	// Обнуление внедиагональных элементов нижнего треугльника 
	for (int i = ig[ind]; i < ig[ind + 1]; i++) {
		ggl[i] = 0;
	}

	// Обнуление внедиагональных элементов верхнего 
	for (int i = ig[ind+1]; i < ig[n]; i++) {
		if (jg[i] == ind) {
			ggu[i] = 0;
		}
	}
}
