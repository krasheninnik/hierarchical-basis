#include "task.h"
#include <functional>
#include <algorithm>
#include <assert.h>
#include <iostream>
#include "stdafx.h"
#include "pardiso.h"
#include "FormatConverter.h"

Node::Node(double _x, double _y) {
	x = _x; y = _y;
}

FiniteElem::FiniteElem(std::vector<int>& _nodes) {
	nodes = _nodes;
	info = _nodes; // first enumerate functions assotiated with nodes
}

void FiniteElem::expandInfo() {
	// 9 functions - for second order
	info.resize(9);
	// fill not used now functions to -1 index
	std::fill(info.begin() + 4, info.end(), -1);
}

void Task::init() {
	gaussIntegration.init(gaussIntegrationOrder);

	initSpaceGrid();
	formatingGlobalMatrixPortrait();
	initParams();

	localMatrix = std::vector<std::vector<double>>(numOfBasisFunctions, std::vector<double>(numOfBasisFunctions));
	localRightPart = std::vector<double>(numOfBasisFunctions);

	// ��������� ������ ��� ������� ����
	const int elemsInMatrix = globalMatrix.getAmountElems();
	const int dimMatrix = globalMatrix.getDim();

	f = std::vector<double>(dimMatrix);
	x = std::vector<double>(dimMatrix);
	temp = std::vector<double>(dimMatrix);
	temp2 = std::vector<double>(dimMatrix);

	// init local functions:
	// Hierarchical basis functions
	func1 f1 = [](double x) { return (1 - x) / 2; };
	func1 f2 = [](double x) { return (1 + x) / 2; };
	func1 f3 = [](double x) { return 1 - x * x; };

	// Hierarchical basis derivatives
	func1 df1 = [](double x) { return - 0.5; };
	func1 df2 = [](double x) { return	0.5; };
	func1 df3 = [](double x) { return - 2 * x; };

	localFunc.resize(numOfBasisFunctions);
	localDx.resize(numOfBasisFunctions);
	localDy.resize(numOfBasisFunctions);
	localFuncTemp.resize(numOfBasisFunctions);
	localDxTemp.resize(numOfBasisFunctions);
	localDyTemp.resize(numOfBasisFunctions);

	// 0 1 2 3 8 4 6 5 7 

	localFuncTemp[0] = [f1, f2, f3](double x, double y, double dx, double dy) { return f1(x/dx) * f1(y/dy); };
	localFuncTemp[1] = [f1, f2, f3](double x, double y, double dx, double dy) { return f2(x/dx) * f1(y/dy); };
	localFuncTemp[2] = [f1, f2, f3](double x, double y, double dx, double dy) { return f1(x/dx) * f2(y/dy); };
	localFuncTemp[3] = [f1, f2, f3](double x, double y, double dx, double dy) { return f2(x/dx) * f2(y/dy); };
	localFuncTemp[4] = [f1, f2, f3](double x, double y, double dx, double dy) { return f1(x/dx) * f3(y/dy); };
	localFuncTemp[5] = [f1, f2, f3](double x, double y, double dx, double dy) { return f2(x/dx) * f3(y/dy); };
	localFuncTemp[6] = [f1, f2, f3](double x, double y, double dx, double dy) { return f3(x/dx) * f1(y/dy); };
	localFuncTemp[7] = [f1, f2, f3](double x, double y, double dx, double dy) { return f3(x/dx) * f2(y/dy); };
	localFuncTemp[8] = [f1, f2, f3](double x, double y, double dx, double dy) { return f3(x/dx) * f3(y/dy); };

	localDxTemp[0] = [f1, f2, f3, df1, df2, df3](double x, double y, double dx, double dy) { return df1(x / dx) * f1(y / dy) / dx; };
	localDxTemp[1] = [f1, f2, f3, df1, df2, df3](double x, double y, double dx, double dy) { return df2(x / dx) * f1(y / dy) / dx; };
	localDxTemp[2] = [f1, f2, f3, df1, df2, df3](double x, double y, double dx, double dy) { return df1(x / dx) * f2(y / dy) / dx; };
	localDxTemp[3] = [f1, f2, f3, df1, df2, df3](double x, double y, double dx, double dy) { return df2(x / dx) * f2(y / dy) / dx; };
	localDxTemp[4] = [f1, f2, f3, df1, df2, df3](double x, double y, double dx, double dy) { return df1(x / dx) * f3(y / dy) / dx; };
	localDxTemp[5] = [f1, f2, f3, df1, df2, df3](double x, double y, double dx, double dy) { return df2(x / dx) * f3(y / dy) / dx; };
	localDxTemp[6] = [f1, f2, f3, df1, df2, df3](double x, double y, double dx, double dy) { return df3(x / dx) * f1(y / dy) / dx; };
	localDxTemp[7] = [f1, f2, f3, df1, df2, df3](double x, double y, double dx, double dy) { return df3(x / dx) * f2(y / dy) / dx; };
	localDxTemp[8] = [f1, f2, f3, df1, df2, df3](double x, double y, double dx, double dy) { return df3(x / dx) * f3(y / dy) / dx; };

	localDyTemp[0] = [f1, f2, f3, df1, df2, df3](double x, double y, double dx, double dy) { return f1(x / dx) * df1(y / dy) / dy; };
	localDyTemp[1] = [f1, f2, f3, df1, df2, df3](double x, double y, double dx, double dy) { return f2(x / dx) * df1(y / dy) / dy; };
	localDyTemp[2] = [f1, f2, f3, df1, df2, df3](double x, double y, double dx, double dy) { return f1(x / dx) * df2(y / dy) / dy; };
	localDyTemp[3] = [f1, f2, f3, df1, df2, df3](double x, double y, double dx, double dy) { return f2(x / dx) * df2(y / dy) / dy; };
	localDyTemp[4] = [f1, f2, f3, df1, df2, df3](double x, double y, double dx, double dy) { return f1(x / dx) * df3(y / dy) / dy; };
	localDyTemp[5] = [f1, f2, f3, df1, df2, df3](double x, double y, double dx, double dy) { return f2(x / dx) * df3(y / dy) / dy; };
	localDyTemp[6] = [f1, f2, f3, df1, df2, df3](double x, double y, double dx, double dy) { return f3(x / dx) * df1(y / dy) / dy; };
	localDyTemp[7] = [f1, f2, f3, df1, df2, df3](double x, double y, double dx, double dy) { return f3(x / dx) * df2(y / dy) / dy; };
	localDyTemp[8] = [f1, f2, f3, df1, df2, df3](double x, double y, double dx, double dy) { return f3(x / dx) * df3(y / dy) / dy; };

	///////////
	localFunc[0] = localFuncTemp[0];
	localFunc[1] = localFuncTemp[1];
	localFunc[2] = localFuncTemp[2];
	localFunc[3] = localFuncTemp[3];
	localFunc[4] = localFuncTemp[8];
	localFunc[5] = localFuncTemp[4];
	localFunc[6] = localFuncTemp[6];
	localFunc[7] = localFuncTemp[5];
	localFunc[8] = localFuncTemp[7];

	localDx[0] = localDxTemp[0];
	localDx[1] = localDxTemp[1];
	localDx[2] = localDxTemp[2];
	localDx[3] = localDxTemp[3];
	localDx[4] = localDxTemp[8];
	localDx[5] = localDxTemp[4];
	localDx[6] = localDxTemp[6];
	localDx[7] = localDxTemp[5];
	localDx[8] = localDxTemp[7];
									   
	localDy[0] = localDyTemp[0];
	localDy[1] = localDyTemp[1];
	localDy[2] = localDyTemp[2];
	localDy[3] = localDyTemp[3];
	localDy[4] = localDyTemp[8];
	localDy[5] = localDyTemp[4];
	localDy[6] = localDyTemp[6];
	localDy[7] = localDyTemp[5];
	localDy[8] = localDyTemp[7];
	//////////
}

void Task::fillAxisGrid(std::vector<double>& axis, double a, double b, int steps, double coef, const int k) {
	double step = 0;
	double point = a;

	if (abs(coef - 1) > 1e-13) step = (b - a) * (1 - coef) / (1 - pow(coef, steps));
	else step = (b - a) / steps;

	// ���������� ������� ����
	coef = pow(coef, 1.0 / k);
	steps *= k;

	double stepsCoef = 0;
	for (int i = 0; i < k; i++) stepsCoef += pow(coef, i);
	step /= stepsCoef;
	
	if(axis.empty()) axis.push_back(point);

	for (int i = 0; i < steps; i++) {
		point += step;
		step *= coef;
		axis.push_back(point);
	}
}

void Task::initSpaceGrid() {
	std::fstream fin(R"(input\grid_space.txt)");
	gridDiv = 1;
	fin >> gridDiv;
	int k = pow(2, gridDiv - 1);

	double startPoint, numOfSegm, step;
	std::vector<int> numOfElems;
	double a, b, steps, coef;

	std::vector<double> xAxisGrid;
	std::vector<double> yAxisGrid;

	fin >> numOfSegm;
	for (int i = 0; i < numOfSegm; i++) {
		fin >> a >> b >> steps >> coef;
		fillAxisGrid(xAxisGrid, a, b, steps, coef, k);
		fillAxisGrid(xAxisWithoutDiv, a, b, steps, coef, 1);
	}

	fin >> numOfSegm;
	for (int i = 0; i < numOfSegm; i++) {
		fin >> a >> b >> steps >> coef;
		fillAxisGrid(yAxisGrid, a, b, steps, coef, k);
		fillAxisGrid(yAxisWithoutDiv, a, b, steps, coef, 1);
	}

	for (double y : yAxisGrid) {
		for (double x : xAxisGrid) {
			nodes.push_back(Node(x, y));
		}
	}

	xaxis = xAxisGrid;
	yaxis = yAxisGrid;

	// ������������ ������� ���������
	const int NODE_POINTS_SIZE = 4;
	std::vector<int> globalNodes = std::vector<int>(NODE_POINTS_SIZE);

	// ���������� ��������� � * ��������� � ����� z-������
	nx = xAxisGrid.size() - 1;
	ny = yAxisGrid.size() - 1;

	// ���������� ����� � * ��������� � ����� z-������
	npx = xAxisGrid.size();
	npy = yAxisGrid.size();

	// ���������� ������� ���������
	int firstElemNode = 0;
	for (int yi = 0; yi < yAxisGrid.size() - 1; yi++) {
		for (int xi = 0; xi < xAxisGrid.size() - 1; xi++) {
			// it's just XY points to determine element
			globalNodes[0] = firstElemNode;
			globalNodes[1] = firstElemNode + 1;
			globalNodes[2] = firstElemNode + npx;
			globalNodes[3] = firstElemNode + 1 + npx;

			firstElemNode++;
			elems.push_back(FiniteElem(globalNodes));
		}
		firstElemNode++;
	}


	// second order elems:
	
	auto getElemInd = [this](int xi, int yi) {return yi * this->nx + xi; };
	auto getLeftNeighInd = [this](int xi, int yi) {return yi * this->nx + xi - 1; };
	auto getRightNeighInd = [this](int xi, int yi) {return yi * this->nx + xi + 1; };
	auto getUpNeighInd = [this](int xi, int yi) {return (yi + 1) * this->nx + xi; };
	auto getDownNeighInd = [this](int xi, int yi) {return (yi - 1) * this->nx + xi; };

	std::vector<std::tuple<int, int, int, int>> secondOrderElemsAreas;
	int secondOrderElemsAreasSize = 0;
	fin >> secondOrderElemsAreasSize;
	secondOrderElemsAreas.resize(secondOrderElemsAreasSize);
	for (int i = 0; i < secondOrderElemsAreasSize; i++) {
		int x0 = 0, y0 = 0, x1 = 0, y1 = 0;
		fin >> x0 >> y0 >> x1 >> y1;

		// consider grid division:
		x0 *= k;
		y0 *= k;
		x1 = (x1 + 1) * k - 1;
		y1 = (y1 + 1) * k - 1;
		secondOrderElemsAreas[i] = std::make_tuple(x0, y0, x1, y1);

		// mark elems with second order:
		for (int yi = y0; yi <= y1; yi++) {
			for (int xi = x0; xi <= x1; xi++) {
				elems[getElemInd(xi, yi)].functionOrder = 2;
				elems[getElemInd(xi, yi)].expandInfo();
			}
		}
	}

	// set lambda and gamma:
	int areasSize = 0;
	double theLambda = 0, theGamma = 0;
	fin >> areasSize;
	for (int i = 0; i < areasSize; i++) {
		int x0 = 0, y0 = 0, x1 = 0, y1 = 0;
		fin >> x0 >> y0 >> x1 >> y1 >> theLambda >> theGamma;

		// consider grid division:
		x0 *= k;
		y0 *= k;
		x1 = (x1 + 1) * k - 1;
		y1 = (y1 + 1) * k - 1;

		//  set lambda and gamma
		for (int yi = y0; yi <= y1; yi++) {
			for (int xi = x0; xi <= x1; xi++) {
				elems[getElemInd(xi, yi)].lambda = theLambda;
				elems[getElemInd(xi, yi)].gamma = theGamma;
			}
		}
	}

	// check setting lambda and gamma:
	for (auto& el : elems) {
		assert(el.lambda != 0 || el.gamma != 0);
	}
	// fill "info" array for second elems 
	int lastFunction = nodes.size() - 1;
	auto setElemInfoIfItNotInited = [this, &lastFunction](int sourceElemInd, int sourceInfoInd) {
		if (elems[sourceElemInd].info[sourceInfoInd] == -1) elems[sourceElemInd].info[sourceInfoInd] = ++lastFunction;
	};

	auto setElemInfoIfItNotInitedWithNeigh = [this, &lastFunction](int sourceElemInd,  int sourceInfoInd, int neighElemInd, int neighInfoInd) {
		auto& elem = this->elems[sourceElemInd];
		auto& neighElem = this->elems[neighElemInd];
		if (elem.info[sourceInfoInd] == -1 && this->elems[neighElemInd].functionOrder == 2) {
			if (neighElem.info[neighInfoInd] == -1) {
				elem.info[sourceInfoInd] = ++lastFunction;
				neighElem.info[neighInfoInd] = elem.info[sourceInfoInd];
			}
			else {
				elem.info[sourceInfoInd] = neighElem.info[neighInfoInd];
			}
		}
	};

	for (int i = 0; i < secondOrderElemsAreas.size(); i++) {
		auto area = secondOrderElemsAreas[i];
		int x0 = std::get<0>(area);
		int y0 = std::get<1>(area);
		int x1 = std::get<2>(area);
		int y1 = std::get<3>(area);


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


		for (int yi = y0; yi <= y1; yi++) {
			// if only one lvl by y:
			if (yi == 0 && yi == ny - 1) {
				// there only one lvl by x:
				if (nx - 1 == 0) {
					for (int i = 4; i < 9; i++) {
						elems[0].info[i] = ++lastFunction;
					}
					continue;
				}

				for (int xi = x0; xi <= x1; xi++) {
					// check left elems:
					if (xi == 0) {
						int sourceElemInd = getElemInd(xi, yi);
						int rightNeighElemInd = getRightNeighInd(xi, yi);

						setElemInfoIfItNotInited(sourceElemInd, 4);
						setElemInfoIfItNotInited(sourceElemInd, 5);
						setElemInfoIfItNotInited(sourceElemInd, 6);
						setElemInfoIfItNotInitedWithNeigh(sourceElemInd, 7, rightNeighElemInd, 5);
						setElemInfoIfItNotInited(sourceElemInd, 8);
						continue;
					}

					// check right elems
					if (xi == nx - 1) {
						int sourceElemInd = getElemInd(xi, yi);
						int leftNeighElemInd = getLeftNeighInd(xi, yi);

						setElemInfoIfItNotInited(sourceElemInd, 4);
						setElemInfoIfItNotInitedWithNeigh(sourceElemInd, 5, leftNeighElemInd, 7);
						setElemInfoIfItNotInited(sourceElemInd, 6);
						setElemInfoIfItNotInited(sourceElemInd, 7);
						setElemInfoIfItNotInited(sourceElemInd, 8);
						continue;
					}

					// if elems in the middle:
					int sourceElemInd = getElemInd(xi, yi);
					int leftNeighElemInd = getLeftNeighInd(xi, yi);
					int rightNeighElemInd = getRightNeighInd(xi, yi);

					setElemInfoIfItNotInited(sourceElemInd, 4);
					setElemInfoIfItNotInitedWithNeigh(sourceElemInd, 5, leftNeighElemInd, 7);
					setElemInfoIfItNotInited(sourceElemInd, 6);
					setElemInfoIfItNotInitedWithNeigh(sourceElemInd, 7, rightNeighElemInd, 5);
					setElemInfoIfItNotInited(sourceElemInd, 8);
				}
				continue;
			}

			// check bottom elems:
			if (yi == 0) {
				// there only one lvl by x:
				if (nx - 1 == 0) {
					int sourceElemInd = 0;
					int upNeighElemInd = 1;

					setElemInfoIfItNotInited(sourceElemInd, 4);
					setElemInfoIfItNotInited(sourceElemInd, 5);
					setElemInfoIfItNotInited(sourceElemInd, 6);
					setElemInfoIfItNotInited(sourceElemInd, 7);
					setElemInfoIfItNotInitedWithNeigh(sourceElemInd, 8, upNeighElemInd, 6);

					continue;
				}
				else {
					for (int xi = x0; xi <= x1; xi++) {
						// check left elems:
						if (xi == 0) {
							int sourceElemInd = getElemInd(xi, yi);
							int rightNeighElemInd = getRightNeighInd(xi, yi);
							int upNeighElemInd = getUpNeighInd(xi, yi);

							setElemInfoIfItNotInited(sourceElemInd, 4);
							setElemInfoIfItNotInited(sourceElemInd, 5);
							setElemInfoIfItNotInited(sourceElemInd, 6);
							setElemInfoIfItNotInitedWithNeigh(sourceElemInd, 7, rightNeighElemInd, 5);
							setElemInfoIfItNotInitedWithNeigh(sourceElemInd, 8, upNeighElemInd, 6);
							continue;
						}

						// check right elems
						if (xi == nx - 1) {
							int sourceElemInd = getElemInd(xi, yi);
							int leftNeighElemInd = getLeftNeighInd(xi, yi);
							int upNeighElemInd = getUpNeighInd(xi, yi);

							setElemInfoIfItNotInited(sourceElemInd, 4);
							setElemInfoIfItNotInitedWithNeigh(sourceElemInd, 5, leftNeighElemInd, 7);
							setElemInfoIfItNotInited(sourceElemInd, 6);
							setElemInfoIfItNotInited(sourceElemInd, 7);
							setElemInfoIfItNotInitedWithNeigh(sourceElemInd, 8, upNeighElemInd, 6);
							continue;
						}

						// if elems in the middle:
						int sourceElemInd = getElemInd(xi, yi);
						int rightNeighElemInd = getRightNeighInd(xi, yi);
						int leftNeighElemInd = getLeftNeighInd(xi, yi);
						int upNeighElemInd = getUpNeighInd(xi, yi);

						setElemInfoIfItNotInited(sourceElemInd, 4);
						setElemInfoIfItNotInitedWithNeigh(sourceElemInd, 5, leftNeighElemInd, 7);
						setElemInfoIfItNotInited(sourceElemInd, 6);
						setElemInfoIfItNotInitedWithNeigh(sourceElemInd, 7, rightNeighElemInd, 5);
						setElemInfoIfItNotInitedWithNeigh(sourceElemInd, 8, upNeighElemInd, 6);
						continue;
					}
				}
				continue;
			}

			// check top elems
			if (yi == ny - 1) {
				// there only one lvl by x:
				if (nx - 1 == 0) {
					int sourceElemInd = yi;
					int downNeighElemInd = yi - 1;

					setElemInfoIfItNotInited(sourceElemInd, 4);
					setElemInfoIfItNotInited(sourceElemInd, 5);
					setElemInfoIfItNotInitedWithNeigh(sourceElemInd, 6, downNeighElemInd, 8);
					setElemInfoIfItNotInited(sourceElemInd, 7);
					setElemInfoIfItNotInited(sourceElemInd, 8);
					continue;
				}
				else {
					for (int xi = x0; xi <= x1; xi++) {
						// check left elems:
						if (xi == 0) {
							int sourceElemInd = getElemInd(xi, yi);
							int rightNeighElemInd = getRightNeighInd(xi, yi);
							int downNeighElemInd = getDownNeighInd(xi, yi);

							setElemInfoIfItNotInited(sourceElemInd, 4);
							setElemInfoIfItNotInited(sourceElemInd, 5);
							setElemInfoIfItNotInitedWithNeigh(sourceElemInd, 6, downNeighElemInd, 8);
							setElemInfoIfItNotInitedWithNeigh(sourceElemInd, 7, rightNeighElemInd, 5);
							setElemInfoIfItNotInited(sourceElemInd, 8);
							continue;
						}

						// check right elems
						if (xi == nx - 1) {
							int sourceElemInd = getElemInd(xi, yi);
							int leftNeighElemInd = getLeftNeighInd(xi, yi);
							int downNeighElemInd = getDownNeighInd(xi, yi);

							setElemInfoIfItNotInited(sourceElemInd, 4);
							setElemInfoIfItNotInitedWithNeigh(sourceElemInd, 5, leftNeighElemInd, 7);
							setElemInfoIfItNotInitedWithNeigh(sourceElemInd, 6, downNeighElemInd, 8);
							setElemInfoIfItNotInited(sourceElemInd, 7);
							setElemInfoIfItNotInited(sourceElemInd, 8);
							continue;
						}

						// if elems in the middle:
						int sourceElemInd = getElemInd(xi, yi);
						int rightNeighElemInd = getRightNeighInd(xi, yi);
						int leftNeighElemInd = getLeftNeighInd(xi, yi);
						int downNeighElemInd = getDownNeighInd(xi, yi);

						setElemInfoIfItNotInited(sourceElemInd, 4);
						setElemInfoIfItNotInitedWithNeigh(sourceElemInd, 5, leftNeighElemInd, 7);
						setElemInfoIfItNotInitedWithNeigh(sourceElemInd, 6, downNeighElemInd, 8);
						setElemInfoIfItNotInitedWithNeigh(sourceElemInd, 7, rightNeighElemInd, 5);
						setElemInfoIfItNotInited(sourceElemInd, 8);
						continue;
					}
					continue;
				}
			}

			// check elems in the middle
			for (int xi = x0; xi <= x1; xi++) {
				// there only one lvl by x:
				if (nx - 1 == 0) {
					int sourceElemInd = yi;
					int upNeighElemInd = yi + 1;
					int downNeighElemInd = yi + -1;

					setElemInfoIfItNotInited(sourceElemInd, 4);
					setElemInfoIfItNotInited(sourceElemInd, 5);
					setElemInfoIfItNotInitedWithNeigh(sourceElemInd, 6, downNeighElemInd, 8);
					setElemInfoIfItNotInited(sourceElemInd, 7);
					setElemInfoIfItNotInitedWithNeigh(sourceElemInd, 8, upNeighElemInd, 6);
					continue;
				}
				else {
					for (int xi = x0; xi <= x1; xi++) {
						// check left elems:
						if (xi == 0) {
							int sourceElemInd = getElemInd(xi, yi);
							int rightNeighElemInd = getRightNeighInd(xi, yi);
							int downNeighElemInd = getDownNeighInd(xi, yi);
							int upNeighElemInd = getUpNeighInd(xi, yi);

							setElemInfoIfItNotInited(sourceElemInd, 4);
							setElemInfoIfItNotInited(sourceElemInd, 5);
							setElemInfoIfItNotInitedWithNeigh(sourceElemInd, 6, downNeighElemInd, 8);
							setElemInfoIfItNotInitedWithNeigh(sourceElemInd, 7, rightNeighElemInd, 5);
							setElemInfoIfItNotInitedWithNeigh(sourceElemInd, 8, upNeighElemInd, 6);
							continue;
						}

						// check right elems
						if (xi == nx - 1) {
							int sourceElemInd = getElemInd(xi, yi);
							int leftNeighElemInd = getLeftNeighInd(xi, yi);
							int downNeighElemInd = getDownNeighInd(xi, yi);
							int upNeighElemInd = getUpNeighInd(xi, yi);

							setElemInfoIfItNotInited(sourceElemInd, 4);
							setElemInfoIfItNotInitedWithNeigh(sourceElemInd, 5, leftNeighElemInd, 7);
							setElemInfoIfItNotInitedWithNeigh(sourceElemInd, 6, downNeighElemInd, 8);
							setElemInfoIfItNotInited(sourceElemInd, 7);
							setElemInfoIfItNotInitedWithNeigh(sourceElemInd, 8, upNeighElemInd, 6);
							continue;
						}

						// if elems in the middle:
						int sourceElemInd = getElemInd(xi, yi);
						int rightNeighElemInd = getRightNeighInd(xi, yi);
						int leftNeighElemInd = getLeftNeighInd(xi, yi);
						int downNeighElemInd = getDownNeighInd(xi, yi);
						int upNeighElemInd = getUpNeighInd(xi, yi);

						setElemInfoIfItNotInited(sourceElemInd, 4);
						setElemInfoIfItNotInitedWithNeigh(sourceElemInd, 5, leftNeighElemInd, 7);
						setElemInfoIfItNotInitedWithNeigh(sourceElemInd, 6, downNeighElemInd, 8);
						setElemInfoIfItNotInitedWithNeigh(sourceElemInd, 7, rightNeighElemInd, 5);
						setElemInfoIfItNotInitedWithNeigh(sourceElemInd, 8, upNeighElemInd, 6);
						continue;
					}
					continue;
				}
			}
		}
	}

	

	// forming boundariesElems vectors
	boundariesElemsTop.resize(nx);
	boundariesElemsBottom.resize(nx);
	boundariesElemsLeft.resize(ny);
	boundariesElemsRight.resize(ny);

	for (int i = 0; i < nx; i++) boundariesElemsBottom[i] = i;
	for (int i = 0, elInd = nx*(ny-1); i < nx; i++, elInd++) boundariesElemsTop[i] = elInd;
	for (int i = 0, elInd = 0; i < ny; i++, elInd+=nx) boundariesElemsLeft[i] = elInd;
	for (int i = 0, elInd = nx - 1; i < ny; i++, elInd += nx) boundariesElemsRight[i] = elInd;

	// for sake of debug elems orders:
	bool drawOrders = false;
	bool drawElems = false;

	if(drawOrders) {

		std::cout << "####|";
		for (int xi = 0; xi < nx; xi++) {
			if (xi >= 10) {
				std::cout << xi << "| ";
			}
			else {
				std::cout << xi << " | ";
			}
		}
		std::cout << std::endl;

		for (int yi = 0; yi < ny; yi++) {
			if (yi < 10) std::cout << " ";
			std::cout << yi << ": |";
			for (int xi = 0; xi < nx; xi++) {
				int order = elems[yi * nx + xi].functionOrder;
				if (order == 1) {
					std::cout << "-" << " | ";

				}
				else {
					std::cout << order << " | ";
				}


			}
			std::cout << std::endl;
		}

	}

	if (drawElems) {
		for (int i = 0; i < secondOrderElemsAreas.size(); i++) {
			auto area = secondOrderElemsAreas[i];
			int x0 = std::get<0>(area);
			int y0 = std::get<1>(area);
			int x1 = std::get<2>(area);
			int y1 = std::get<3>(area);
			
			std::cout << "x0,y0,x1,y1: " << x0 << " " << y0 << " " << x1 << " " << y1 << std::endl;

			for (int yi = y0; yi <= y1; yi++) {

				for (int xi = x0; xi <= x1; xi++) {
					auto& elem = elems[getElemInd(xi, yi)];
					std::cout << "[" << setw(3) << elem.info[0] << "] [" << setw(3) << elem.info[6] << "] [" << setw(3) << elem.info[1] << "] |";
				}
				std::cout << std::endl;

				for (int xi = x0; xi <= x1; xi++) {
					auto& elem = elems[getElemInd(xi, yi)];
					std::cout << "[" << setw(3) << elem.info[5] << "] [" << setw(3) << elem.info[4] << "] [" << setw(3) << elem.info[7] << "] |";
				}
				std::cout << std::endl;

				for (int xi = x0; xi <= x1; xi++) {
					auto& elem = elems[getElemInd(xi, yi)];
					std::cout << "[" << setw(3) <<  elem.info[2] << "] [" << setw(3) << elem.info[8] << "] [" << setw(3) << elem.info[3] << "] |";
				}
		
				std::cout << std::endl << "----------" << std::endl;

			}
			std::cout << std::endl;
		}
	}


	fin.close();
}
	
void Task::formatingGlobalMatrixPortrait() {
	//DIM = elems.back().info.back() + 1;
	for (auto el : elems) {
		int maximum = *max_element(el.info.begin(), el.info.end());
		if (maximum > DIM) DIM = maximum;
	}
	DIM++;
	std::cout << "Matrix DIM: " << DIM << std::endl;


	std::vector<std::vector<int>> temp(DIM);

	std::function<bool(std::vector<int>&, double)> hasValue = [](std::vector<int>& v, double value) {
		if (std::find(v.begin(), v.end(), value) == v.end()) return false;
		else return true;
	};

	// ������������ ������
	for (auto& el : elems) {
		/*
			[2] - [8] - [3]
			 |	   |	 |
			[5] - [4] - [7]
			 |	   |	 |
			[0] - [6] - [1]
		*/
		
		auto elInfo = el.info;

		// delete "-1 functions (fictious)"
		elInfo.erase(std::remove_if(elInfo.begin(), elInfo.end(), [](int x) {return x == -1; }), elInfo.end());
		
		// increase sort functions
		std::sort(elInfo.begin(), elInfo.end());

		for (int rowInd = 1; rowInd < elInfo.size(); rowInd++) {
			for (int columnInd = 0; columnInd < rowInd; columnInd++) {
				if (!hasValue(temp[elInfo[rowInd]], elInfo[columnInd])) temp[elInfo[rowInd]].push_back(elInfo[columnInd]);
			}
		}
	}

	// ���������� ��������
	for (auto& columns : temp) {
		std::sort(columns.begin(), columns.end());
	}

	// ���������� jj �������:
	std::vector<int> jg;
	for (auto& vec : temp) {
		std::copy(vec.begin(), vec.end(), back_inserter(jg));
	}

	// ���������� ig �������
	std::vector<int> ig;
	ig.push_back(0);
	for (auto& vec : temp) {
		ig.push_back(ig.back() + vec.size());
	}

	// ������������� ���������� �������
	globalMatrix.init(DIM, ig, jg);
}

void Task::calculateLocalMatrix(const FiniteElem &el) {
	double x1 = nodes[el.nodes[1]].x;
	double x0 = nodes[el.nodes[0]].x;

	double y1 = nodes[el.nodes[2]].y;
	double y0 = nodes[el.nodes[0]].y;

	double dx = x1 - x0;
	double dy = y1 - y0;

	assert(dx > 0);
	assert(dy > 0);

	const double lambda = el.lambda;
	const double gamma = el.gamma;

	// calculate (gradU gradV) integral
	for (int row = 0; row < el.info.size(); row++) {
		for (int column = 0; column <= row; column++) {
			auto func = [this, row, column, dx, dy, lambda, gamma](double x, double y) {return lambda * this->localDx[row](x, y, dx, dy) * this->localDx[column](x, y, dx, dy) + this->localDy[row](x, y, dx, dy) * this->localDy[column](x, y, dx, dy); };
			double integral = gaussIntegration.nPointsGauss(-dx, dx, -dy, dy, func) * 4;
			localMatrix[row][column] = integral;
		}
	}

	// calculate (U V) integral 
	for (int row = 0; row < el.info.size(); row++) {
		for (int column = 0; column <= row; column++) {
			auto func = [this, row, column, dx, dy, lambda, gamma](double x, double y) {return gamma * this->localFunc[row](x, y, dx, dy) * this->localFunc[column](x, y, dx, dy); };
			double integral = gaussIntegration.nPointsGauss(-dx, dx, -dy, dy, func);
			localMatrix[row][column] += integral;
		}
	}

	int c = 3;
}

void Task::calculateLocalRightPart(const FiniteElem& el) {
	double x1 = nodes[el.nodes[1]].x;
	double x0 = nodes[el.nodes[0]].x;

	double y1 = nodes[el.nodes[2]].y;
	double y0 = nodes[el.nodes[0]].y;

	double dx = x1 - x0;
	double dy = y1 - y0;

	const double lambda = el.lambda;
	const double gamma = el.gamma;

	// calculate (f V) integral
	for (int row = 0; row < localRightPart.size(); row++) {
		auto func = [this, row, dx, dy, x0, y0, lambda, gamma](double x, double y) {return this->rightPartFunction(x0 + 0.5*(dx + x), y0 + 0.5 * (dy + y), lambda, gamma) * this->localFunc[row](x, y, dx, dy); };
		localRightPart[row] = gaussIntegration.nPointsGauss(-dx, dx, -dy, dy, func);
	}

	int c = 3;
}

void Task::addLocalMatrixToGlobal(const FiniteElem& elem) {
	// add diagonal elements:
	for (int i = 0; i < elem.info.size(); i++) {
		if (elem.info[i] != -1) {
			globalMatrix.addDiagElem(elem.info[i], localMatrix[i][i]);
		}
	}

	// add function links
	for (int i = 1; i < elem.info.size(); i++) {
		if (elem.info[i] != -1) {
			for (int j = 0; j < i; j++)
				if (elem.info[j] != -1) {
					// probably something will go wrong !!!
					int rowGlobalInd = max(elem.info[i], elem.info[j]);
					int columnGlobalInd = min(elem.info[i], elem.info[j]);

					//std::cout << "i: " << setw(4) << rowGlobalInd << " j: " << setw(4) << columnGlobalInd << std::endl;

					globalMatrix.addElem(rowGlobalInd, columnGlobalInd, localMatrix[i][j]);
				}
		}
	}
	//std::cout << "-------------------------" << std::endl;
}

void Task::addLocalRigtPartToGlobal(const FiniteElem& elem) {
	for (int i = 0; i < elem.info.size(); i++) {
		if (elem.info[i] != -1) {
			f[elem.info[i]] += localRightPart[i];
		}
	}
}

void Task::setFirstBoundaryConditions() {
	// =============================================================================================
	// ADD BOTTOM ELEMS:
	//std::cout << "Bottom: " << std::endl;
	for (int ind : boundariesElemsBottom) {
		auto& elem = elems[ind];

		double x1 = nodes[elem.nodes[1]].x;
		double x0 = nodes[elem.nodes[0]].x;

		double y1 = nodes[elem.nodes[2]].y;
		double y0 = nodes[elem.nodes[0]].y;

		// add left bottom elem
		globalMatrix.setFirstBoundaryCondition(elem.info[0]);
		f[elem.info[0]] = boundaryFunction(x0, y0);

		//std::cout << elem.info[0] << std::endl;

		if (elem.functionOrder == 2) {
			// add mid bottom elem
			globalMatrix.setFirstBoundaryCondition(elem.info[6]);
			f[elem.info[6]] = boundaryFunction((x0 + x1) / 2, y0) - 0.5*(boundaryFunction(x0,y0) + boundaryFunction(x1,y0));
			//std::cout << elem.info[6] << std::endl;
		}
	}

	// add last elem
	auto lastElem = elems[boundariesElemsBottom.back()];
	globalMatrix.setFirstBoundaryCondition(lastElem.info[1]);
	f[lastElem.info[1]] = boundaryFunction(nodes[lastElem.nodes[1]].x, nodes[lastElem.nodes[1]].y);
	//std::cout << lastElem.info[1] << std::endl;

	// =============================================================================================
	// ADD TOP ELEMS:
	//std::cout << "Top: " << std::endl;
	for (int ind : boundariesElemsTop) {
		auto& elem = elems[ind];

		double x1 = nodes[elem.nodes[1]].x;
		double x0 = nodes[elem.nodes[0]].x;

		double y1 = nodes[elem.nodes[2]].y;
		double y0 = nodes[elem.nodes[0]].y;

		// add left top elem
		globalMatrix.setFirstBoundaryCondition(elem.info[2]);
		f[elem.info[2]] = boundaryFunction(x0, y1);
		//std::cout << elem.info[2] << std::endl;

		if (elem.functionOrder == 2) {
			// add mid bottom elem
			globalMatrix.setFirstBoundaryCondition(elem.info[8]);
			f[elem.info[8]] = boundaryFunction((x0 + x1) / 2, y1) - 0.5 * (boundaryFunction(x0, y1) + boundaryFunction(x1, y1));
			//std::cout << elem.info[8] << std::endl;
		}
	}

	// add last elem
	lastElem = elems[boundariesElemsTop.back()];
	globalMatrix.setFirstBoundaryCondition(lastElem.info[3]);
	f[lastElem.info[3]] = boundaryFunction(nodes[lastElem.nodes[3]].x, nodes[lastElem.nodes[3]].y);
	//std::cout << lastElem.info[3] << std::endl;

	// =============================================================================================
	// ADD LEFT ELEMS:
	// std::cout << "Left: " << std::endl;
	for (int ind : boundariesElemsLeft) {
		auto& elem = elems[ind];

		double x1 = nodes[elem.nodes[1]].x;
		double x0 = nodes[elem.nodes[0]].x;

		double y1 = nodes[elem.nodes[2]].y;
		double y0 = nodes[elem.nodes[0]].y;

		// add bottom left elem
		globalMatrix.setFirstBoundaryCondition(elem.info[0]);
		f[elem.info[0]] = boundaryFunction(x0, y0);
		//std::cout << elem.info[0] << std::endl;

		if (elem.functionOrder == 2) {
			// add mid bottom elem
			globalMatrix.setFirstBoundaryCondition(elem.info[5]);
			f[elem.info[5]] = boundaryFunction(x0, (y0+y1)/2) - 0.5 * (boundaryFunction(x0, y0) + boundaryFunction(x0, y1));
			//std::cout << elem.info[5] << std::endl;
		}
	}

	// add last elem
	lastElem = elems[boundariesElemsLeft.back()];
	globalMatrix.setFirstBoundaryCondition(lastElem.info[2]);
	f[lastElem.info[2]] = boundaryFunction(nodes[lastElem.nodes[2]].x, nodes[lastElem.nodes[2]].y);
	//std::cout << lastElem.info[2] << std::endl;

	// =============================================================================================
	// ADD RIGHT ELEMS:
	//std::cout << "Right: " << std::endl;
	for (int ind : boundariesElemsRight) {
		auto& elem = elems[ind];

		double x1 = nodes[elem.nodes[1]].x;
		double x0 = nodes[elem.nodes[0]].x;

		double y1 = nodes[elem.nodes[2]].y;
		double y0 = nodes[elem.nodes[0]].y;

		// add bottom right elem
		globalMatrix.setFirstBoundaryCondition(elem.info[1]);
		f[elem.info[1]] = boundaryFunction(x1, y0);
		//std::cout << elem.info[1] << std::endl;

		if (elem.functionOrder == 2) {
			// add mid bottom elem
			globalMatrix.setFirstBoundaryCondition(elem.info[7]);
			f[elem.info[7]] = boundaryFunction(x1, (y0 + y1) / 2) - 0.5 * (boundaryFunction(x1, y0) + boundaryFunction(x1, y1));;
			//std::cout << elem.info[7] << std::endl;
		}
	}

	// add last elem
	lastElem = elems[boundariesElemsRight.back()];
	globalMatrix.setFirstBoundaryCondition(lastElem.info[3]);
	f[lastElem.info[3]] = boundaryFunction(nodes[lastElem.nodes[3]].x, nodes[lastElem.nodes[3]].y);
	//std::cout << lastElem.info[3] << std::endl;
	// =============================================================================================
}

void Task::PARDISOsolve() {
	using namespace std;

	int* ig, * jg, N;
	double* ggu, * ggl, * di, * pr, * _x;

	N = globalMatrix.getDim();
	int koljg = globalMatrix.getAmountElems();

	ig = new int[N + 1];
	jg = new int[koljg];
	di = new double[N];
	ggu = new double[koljg];
	ggl = new double[koljg];

	_x = new double[N];
	pr = new double[N];

	globalMatrix.getIG(ig);
	globalMatrix.getJG(jg);
	globalMatrix.getDI(di);
	globalMatrix.getGGL(ggl);
	globalMatrix.getGGU(ggu);
	for (int i = 0; i < x.size(); i++) _x[i] = x[i];
	for (int i = 0; i < f.size(); i++) pr[i] = f[i];

	cout << "SolveSlae" << endl;

	pardiso_solver prds;

	prds.factorize(N, ig, jg, ggl, ggu, di, 1);

	prds.solve_nrhs(1, pr, _x);
	prds.stop_solver();
	prds.~pardiso_solver();

	for (int i = 0; i < x.size(); i++) x[i] = _x[i];

	delete[] ig ;
	delete[] jg ;
	delete[] di ;
	delete[] ggu;
	delete[] ggl;

	delete[] _x ;
	delete[] pr ;

	return;
}

double Task::resultInXY(double xCord, double yCord) {
	if (xCord < nodes[elems[boundariesElemsBottom[0]].nodes[0]].x ||
		xCord > nodes[elems[boundariesElemsBottom.back()].nodes[1]].x) {
		std::cout << "x = " << xCord << " not in defined area" << std::endl;
		return -999999;
	}

	if (yCord < nodes[elems[boundariesElemsLeft[0]].nodes[0]].y ||
		yCord > nodes[elems[boundariesElemsLeft.back()].nodes[2]].y) {
		std::cout << "y = " << yCord << " not in defined area" << std::endl;
		return -999999;
	}


	int xi = 0;
	for (; xi < boundariesElemsBottom.size(); xi++) {
		auto& el = elems[boundariesElemsBottom[xi]];
		if (nodes[el.nodes[0]].x <= xCord && xCord <= nodes[el.nodes[1]].x) {
			break;
		}
	}

	int yi = 0;
	for (; yi < boundariesElemsLeft.size(); yi++) {
		auto& el = elems[boundariesElemsLeft[yi]];
		if (nodes[el.nodes[0]].y <= yCord && yCord <= nodes[el.nodes[2]].y) {
			break;
		}
	}

	int elemInd = yi * nx + xi;
	auto& elem = elems[elemInd];
	

	double x1 = nodes[elem.nodes[1]].x;
	double x0 = nodes[elem.nodes[0]].x;

	double y1 = nodes[elem.nodes[2]].y;
	double y0 = nodes[elem.nodes[0]].y;

	double dx = x1 - x0;
	double dy = y1 - y0;

	auto changeVariable = [](double x, double a, double b) {
		double halflenght = (b - a) / 2;
		return (x - a - halflenght) / halflenght;
	};

	double res = 0;

	for(int ind = 0; ind < elem.info.size(); ind++) {
		if (elem.info[ind] != -1) {
			double w = x[elem.info[ind]];
			double lf = localFunc[ind](
				changeVariable(xCord, x0, x1),
				changeVariable(yCord, y0, y1),
				1,
				1);

			res += x[elem.info[ind]] * lf;
		}
	}
	
	return res;
}

void Task::initParams() {
	int CASE = 7;

	switch (CASE) {
	case 1: {
		boundaryFunction = [](double x, double y) {return x + y; };
		rightPartFunction = [this](double x, double y, double lambda, double gamma) {return gamma * this->boundaryFunction(x,y); };

		break;
	}
	case 2: {
		boundaryFunction = [](double x, double y) {return  x * x + y * y; };
		rightPartFunction = [this](double x, double y, double lambda, double gamma) {return -4 * lambda + gamma * this->boundaryFunction(x, y); };

		break;
	}
	case 3: {
		boundaryFunction = [](double x, double y) {return ( y * y * y + x * x * x); };
		rightPartFunction = [this](double x, double y, double lambda, double gamma) {return -1 * lambda * (6 * y + 6 * x) + gamma * this->boundaryFunction(x, y); };

		break;
	}
	case 4: {
		boundaryFunction = [](double x, double y) {return (y * y * y * y +  x * x * x * x); };
		rightPartFunction = [this](double x, double y, double lambda, double gamma) {return -1 * (12*y*y + 12 *x * x) * lambda + gamma * this->boundaryFunction(x, y); };

		break;
	}
	case 5: {
		boundaryFunction = [](double x, double y) {return (y * y * y * y * y + x * x * x * x * x); };
		rightPartFunction = [this](double x, double y, double lambda, double gamma) {return -1 * (20 * y * y * y + 20 * x * x * x) * lambda + gamma * this->boundaryFunction(x, y); };

		break;
	}
	case 6: {
		boundaryFunction = [](double x, double y) {return sin(x) * cos(y); };
		rightPartFunction = [this](double x, double y, double lambda, double gamma) {return -1 * (-2* cos(y) *sin(x)) * lambda + gamma * this->boundaryFunction(x, y); };

		break;
	}
	case 7: {
		boundaryFunction = [](double x, double y) {return 3*x + 2*y; };
		rightPartFunction = [this](double x, double y, double lambda, double gamma) {return gamma * this->boundaryFunction(x, y); };

		break;
	}
	default: {
		boundaryFunction = [](double x, double y) {return 10 * x + y; };
		rightPartFunction = [](double x, double y, double lambda, double gamma) {return 10 * x + y; };
		break;
	}
	}


	
};

void Task::solve() {
	for (auto& el : elems) {
		calculateLocalMatrix(el);
		calculateLocalRightPart(el);

		addLocalMatrixToGlobal(el);
		addLocalRigtPartToGlobal(el);
	}

	// ������������ ��������� �� ������� ������������ ������� � �������
	globalMatrix.fillGGU();
	setFirstBoundaryConditions();

	PARDISOsolve();
	
	
	std::cout << setw(4) << "y" << " " << setw(4) << "x" << ":    " << setw(15) << "res" << " " << setw(15) << "exact" << " " << setw(15) << "res - exact" << std::endl;
	std::cout << "--------------------------------------------------------------------------" << std::endl;;

	double xstep = xaxis[1] - xaxis[0];
	double ystep = yaxis[1] - yaxis[0];
	double xMax = xaxis.back();
	double yMax = yaxis.back();

	for (double y : yAxisWithoutDiv) {
		for (double x : xAxisWithoutDiv) {
			double res = resultInXY(x, y);
			double exact = boundaryFunction(x, y);

			std::cout << setw(4) << y << " " << setw(4) << x << ":    " << setw(15) << res << " " << setw(15) << exact << " " << "| abs:" << setw(15) << res - exact << "| ref: " << setw(15) << (res - exact) / (exact) << std::endl;
		}
		std::cout << "--------------------------------------------------------------------------" << std::endl;;
	}

	std::cout << std::endl << "Error norm: " << calculateErrorNorm() << std::endl;
	int b = 1;
}


double Task::calculateErrorNorm() {
	
	double sum = 0;
	/*
	for (double y : yAxisWithoutDiv) {
		for (double x : xAxisWithoutDiv) {
			double res = resultInXY(x, y);
			double exact = boundaryFunction(x, y);

			double error = 0;
			if (exact != 0) {
				error = (exact - res) / (exact);
			}
			sum += error * error;
		}
	}
	double ans = sqrt(sum);
		*/
	double errorSum = 0;
	double exactSum = 0;

	/*
	for (double y : yAxisWithoutDiv) {
		for (double x : xAxisWithoutDiv) {
			double res = resultInXY(x, y);
			double exact = boundaryFunction(x, y);

			errorSum += (exact - res) * (exact - res);
			exactSum += exact * exact;
		}
	}
	*/
	std::vector<double> rx{0.05, 0.1,	0.13,	0.17,	0.25,	0.25,	0.333,	0.37,	0.4,	0.42,	0.483};
	std::vector<double> ry{0.05, 0.35,	0.2,	0.45,	0.35,	0.1,	0.3,	0.44,	0.2,	0.25,	0.495 };

	for (int i = 0; i < rx.size(); i++) {
		double res = resultInXY(rx[i], ry[i]);
		double exact = boundaryFunction(rx[i], ry[i]);

		errorSum += (exact - res) * (exact - res);
		exactSum += exact * exact;
	}


	double ans = sqrt(errorSum) / sqrt(exactSum);

	return ans;
}
