#pragma once

#include <vector>
#include <fstream>
#include "matrix.h"

// ��� ����������� x � y �������� ��������
struct XY {
	int x;
	int y;
};


struct Node {
	Node() = default;
	Node(double _x, double _y);

	double x;
	double y;
};

struct FiniteElem {
	FiniteElem() = default;
	FiniteElem(std::vector<int>& _nodes);
	void expandInfo();

	int functionOrder = 1;
	std::vector<int> nodes;
	std::vector<int> info; 

	/*	local nodes on the elem:
			(indexes in info)
		1-4 for first order
		1-9 for second order

		[2] - [8] - [3]
		 |	   |	 |
		[5] - [4] - [7]
		 |	   |	 |
		[0] - [6] - [1]
	*/

};

class Task {
public:
	void init();

	void solve();
private:
	void PARDISOsolve();
	void initSpaceGrid();

	void formatingGlobalMatrixPortrait();
	void calculateLocalMatrix(const FiniteElem&);
	void calculateLocalRightPart(const FiniteElem&);
	void addLocalMatrixToGlobal(const FiniteElem&);
	void addLocalRigtPartToGlobal(const FiniteElem&);

	void setFirstBoundaryConditions();
	void fillAxisGrid(std::vector<double>& axis, double a, double b , int steps, double coef, const int k);

private:
	static const int elemsInLocalMatrix = 16;
	static const int elemsInLocalRightPart = 7;

	double lambda; // ����������� ��������� (�������� / �������������)

	std::vector<Node> nodes;
	std::vector<FiniteElem> elems;

	std::vector<int> boundariesValue;

	SparseMatrix globalMatrix;
	std::vector<double> f; // ���������� ������ �����
	std::vector<double> x; // �������
	std::vector<double> localMatrix;
	std::vector<double> localRightPart;

	// ���������� � �����:
	int gridDiv;
	int nx; // ��������� � Ox
	int ny;	// ��������� � Oy
	int nz; // ��������� � Oz
	int npx;
	int npy; 
	int npxy; 
	int fxy; 
	int fxz;
	int fyz;
	int faceNumInZlayer;
	int faceNum;

	std::vector<double> temp;
	std::vector<double> temp2;
};