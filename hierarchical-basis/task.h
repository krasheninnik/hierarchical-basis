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
	GaussIntegration gaussIntegration;
	const int gaussIntegrationOrder = 2;

	// Hierarchical basis functions
	func1 ff1 = [](double x) { return (1 - x)/2; };
	func1 ff2 = [](double x) { return (1 + x)/2; };
	func1 ff3 = [](double x) { return 1 - x * x; };

	func2 f1 = [this](double x, double y) { return ff1(x) * ff1(y); };
	func2 f2 = [this](double x, double y) { return ff2(x) * ff1(y); };
	func2 f3 = [this](double x, double y) { return ff1(x) * ff2(y); };
	func2 f4 = [this](double x, double y) { return ff2(x) * ff2(y); };
	func2 f5 = [this](double x, double y) { return ff1(x) * ff3(y); };
	func2 f6 = [this](double x, double y) { return ff2(x) * ff3(y); };
	func2 f7 = [this](double x, double y) { return ff3(x) * ff1(y); };
	func2 f8 = [this](double x, double y) { return ff3(x) * ff2(y); };
	func2 f9 = [this](double x, double y) { return ff3(x) * ff3(y); };

	static const int elemsInLocalMatrix = 16;
	static const int elemsInLocalRightPart = 7;

	double lambda; //  оэффициент уравнени€ (в€зкость / проницаемость)

	std::vector<Node> nodes;
	std::vector<FiniteElem> elems;

	std::vector<int> boundariesValue;

	SparseMatrix globalMatrix;
	std::vector<double> f; // √лобальна€ права€ часть
	std::vector<double> x; // –ешение
	std::vector<double> localMatrix;
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

	std::vector<double> temp;
	std::vector<double> temp2;
};