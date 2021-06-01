#pragma once

#include <vector>
#include <fstream>
#include "matrix.h"
#include "GaussIntegration.h"

// ƒл€ определени€ x и y элемента скважины
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
	using func1 = std::function<double(double)>;
	using func2 = std::function<double(double, double)>;
	using func4 = std::function<double(double, double, double, double)>;

public:
	void init();
	void initParams();

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
	GaussIntegration gaussIntegration;
	const int gaussIntegrationOrder = 5;

	// Hierarchical basis functions and derivates	
	std::vector<func4> localFunc;
	std::vector<func4> localDx;
	std::vector<func4> localDy;
	std::vector<func4> localFuncTemp;
	std::vector<func4> localDxTemp;
	std::vector<func4> localDyTemp;

	static const int numOfBasisFunctions = 9;

	double lambda = 1; //  оэффициент уравнени€ (в€зкость / проницаемость)
	double gamma = 1;

	std::vector<Node> nodes;
	std::vector<FiniteElem> elems;

	std::vector<int> boundariesValue;

	std::vector<int> boundariesElemsTop;
	std::vector<int> boundariesElemsBottom;
	std::vector<int> boundariesElemsLeft;
	std::vector<int> boundariesElemsRight;
	func2 boundaryFunction;
	func2 rightPartFunction;

	SparseMatrix globalMatrix;
	std::vector<double> f; // √лобальна€ права€ часть
	std::vector<double> x; // –ешение
	double resultInXY(double x, double y);

	std::vector<std::vector<double>> localMatrix;
	std::vector<double> localRightPart;

	// »нформаци€ о сетке:
	int gridDiv;
	int nx; // Ёлементов в Ox
	int ny;	// Ёлементов в Oy
	int nz; // Ёлементов в Oz
	int npx;
	int npy; 
	int npxy; 
	int fxy; 
	int fxz;
	int fyz;
	int faceNumInZlayer;
	int faceNum;

	std::vector<double> xaxis;
	std::vector<double> yaxis;

	std::vector<double> temp;
	std::vector<double> temp2;
};