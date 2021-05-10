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
	std::vector<int> ig;		// ������� ��������� - ������ ������
	std::vector<int> jg;		// ������� �������� ����������
	std::vector<double> di;		// ������������ ��������
	std::vector<double> ggu;	// �������� �������� ������������
	std::vector<double> ggl;	// �������� ������� ������������

	int n;						// �����������
};