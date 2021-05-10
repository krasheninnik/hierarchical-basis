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
	std::fstream fin(R"(input\params.txt)");
	fin >> lambda;
	fin.close();

	initSpaceGrid();
	//formatingGlobalMatrixPortrait();

	localMatrix = std::vector<double>(elemsInLocalMatrix);
	localRightPart = std::vector<double>(elemsInLocalRightPart);

	// Аллокация памяти для решения СЛАУ
	const int elemsInMatrix = globalMatrix.getAmountElems();
	const int dimMatrix = globalMatrix.getDim();

	f = std::vector<double>(dimMatrix);
	x = std::vector<double>(dimMatrix);
	temp = std::vector<double>(dimMatrix);
	temp2 = std::vector<double>(dimMatrix);
}

void Task::fillAxisGrid(std::vector<double>& axis, double a, double b, int steps, double coef, const int k) {
	double step = 0;
	double point = a;

	if (abs(coef - 1) > 1e-13) step = (b - a) * (1 - coef) / (1 - pow(coef, steps));
	else step = (b - a) / steps;

	// Вычисление первого шага
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
	}

	fin >> numOfSegm;
	for (int i = 0; i < numOfSegm; i++) {
		fin >> a >> b >> steps >> coef;
		fillAxisGrid(yAxisGrid, a, b, steps, coef, k);
	}

	for (double y : yAxisGrid) {
		for (double x : xAxisGrid) {
			nodes.push_back(Node(x, y));
		}
	}

	// Формирование вектора элементов
	const int NODE_POINTS_SIZE = 4;
	std::vector<int> globalNodes = std::vector<int>(NODE_POINTS_SIZE);

	// Количество элементов в * плоскости в одном z-уровне
	nx = xAxisGrid.size() - 1;
	ny = yAxisGrid.size() - 1;

	// Количество точек в * плоскости в одном z-уровне
	npx = xAxisGrid.size();
	npy = yAxisGrid.size();

	// Заполнение вектора элементов
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
			
			std::cout << "x0,y0,x1,y1" << x0 << " " << y0 << " " << x1 << " " << y1 << std::endl;

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
	const int DIM = elems.back().info.back() + 1;
	std::vector<std::vector<int>> temp(DIM);

	std::function<bool(std::vector<int>&, double)> hasValue = [](std::vector<int>& v, double value) {
		if (std::find(v.begin(), v.end(), value) == v.end()) return false;
		else return true;
	};

	// Формирование связей
	for (auto& el : elems) {
		// По особенностям нумерации: для всех элементов, 
		// Номер грани 0 < номера грани 2
		// Номер грани 3 < номера грани 1
		// Номер грани 4 < номера грани 5
		if (!hasValue(temp[el.info[2]], el.info[0])) temp[el.info[2]].push_back(el.info[0]); // Рассмотрение fxz
		if (!hasValue(temp[el.info[1]], el.info[3])) temp[el.info[1]].push_back(el.info[3]); // Рассмотрение fyz
		if (!hasValue(temp[el.info[5]], el.info[4])) temp[el.info[5]].push_back(el.info[4]); // Рассмотрение fxy

		// for second speeds:
		if (!hasValue(temp[el.info[8]], el.info[6])) temp[el.info[8]].push_back(el.info[6]); // Рассмотрение fxz
		if (!hasValue(temp[el.info[7]], el.info[9])) temp[el.info[7]].push_back(el.info[9]); // Рассмотрение fyz
		if (!hasValue(temp[el.info[11]], el.info[10])) temp[el.info[11]].push_back(el.info[10]); // Рассмотрение fxy

		// Номер любой грани < номера элемента
		if (!hasValue(temp[el.info[12]], el.info[0])) temp[el.info[12]].push_back(el.info[0]);
		if (!hasValue(temp[el.info[12]], el.info[1])) temp[el.info[12]].push_back(el.info[1]);
		if (!hasValue(temp[el.info[12]], el.info[2])) temp[el.info[12]].push_back(el.info[2]);
		if (!hasValue(temp[el.info[12]], el.info[3])) temp[el.info[12]].push_back(el.info[3]);
		if (!hasValue(temp[el.info[12]], el.info[4])) temp[el.info[12]].push_back(el.info[4]);
		if (!hasValue(temp[el.info[12]], el.info[5])) temp[el.info[12]].push_back(el.info[5]);
		if (!hasValue(temp[el.info[12]], el.info[6])) temp[el.info[12]].push_back(el.info[6]);
		if (!hasValue(temp[el.info[12]], el.info[7])) temp[el.info[12]].push_back(el.info[7]);
		if (!hasValue(temp[el.info[12]], el.info[8])) temp[el.info[12]].push_back(el.info[8]);
		if (!hasValue(temp[el.info[12]], el.info[9])) temp[el.info[12]].push_back(el.info[9]);
		if (!hasValue(temp[el.info[12]], el.info[10])) temp[el.info[12]].push_back(el.info[10]);
		if (!hasValue(temp[el.info[12]], el.info[11])) temp[el.info[12]].push_back(el.info[11]);

		// for saturation things:
		if (!hasValue(temp[el.info[13]], el.info[12])) temp[el.info[13]].push_back(el.info[12]);
		if (!hasValue(temp[el.info[14]], el.info[12])) temp[el.info[14]].push_back(el.info[12]);
		if (!hasValue(temp[el.info[14]], el.info[13])) temp[el.info[14]].push_back(el.info[13]);
	}

	// Сортировка столбцов
	for (auto& columns : temp) {
		std::sort(columns.begin(), columns.end());
	}

	// Вычисление jj вектора:
	std::vector<int> jg;
	for (auto& vec : temp) {
		std::copy(vec.begin(), vec.end(), back_inserter(jg));
	}

	// Вычисление ig вектора
	std::vector<int> ig;
	ig.push_back(0);
	for (auto& vec : temp) {
		ig.push_back(ig.back() + vec.size());
	}

	// Инициализация глобальной матрицы
	globalMatrix.init(DIM, ig, jg);
}

void Task::calculateLocalMatrix(const FiniteElem &el) {
	/*
			!!! CHANGE IT !!!
	*/

	double dx = nodes[el.nodes[1]].x - nodes[el.nodes[0]].x;
	double dy = nodes[el.nodes[3]].y - nodes[el.nodes[0]].y;

	assert(dx > 0);
	assert(dy > 0);

	double V = dx * dy;
	const double dt = 1;
	// speed variables
	localMatrix[0] = lambda * V * 1.0 / 3;
	localMatrix[1] = lambda * V * 1.0 / 6;
	localMatrix[2] = lambda * V * 1.0 / 3;
	localMatrix[3] = lambda * V * 1.0 / 3;
	localMatrix[4] = lambda * V * 1.0 / 6;
}

void Task::calculateLocalRightPart(const FiniteElem& el) {
	/*
		!!! CHANGE IT !!!
	*/
	double dx = nodes[el.nodes[1]].x - nodes[el.nodes[0]].x;
	double dy = nodes[el.nodes[3]].y - nodes[el.nodes[0]].y;

	localRightPart[0] = 0; // 
	localRightPart[1] = 0; //
	localRightPart[2] = 0; // если g == 0

}

void Task::addLocalMatrixToGlobal(const FiniteElem& elem) {
	/*
		!!! CHANGE IT !!!
	*/

	// TODO: consider new mapping local to global;

	globalMatrix.addDiagElem(elem.info[0], localMatrix[3]);
	globalMatrix.addDiagElem(elem.info[2], localMatrix[5]);

	globalMatrix.addDiagElem(elem.info[1], localMatrix[2]);
	globalMatrix.addDiagElem(elem.info[3], localMatrix[0]);

	globalMatrix.addDiagElem(elem.info[4], localMatrix[6]);
	globalMatrix.addDiagElem(elem.info[5], localMatrix[8]);

						// строка  	  // столбец	// элемент
	globalMatrix.addElem(elem.info[1], elem.info[3], localMatrix[1]);
	globalMatrix.addElem(elem.info[2], elem.info[0], localMatrix[4]);
	globalMatrix.addElem(elem.info[5], elem.info[4], localMatrix[7]);

	globalMatrix.addElem(elem.info[6], elem.info[3], localMatrix[9]);
	globalMatrix.addElem(elem.info[6], elem.info[1], localMatrix[10]);
	globalMatrix.addElem(elem.info[6], elem.info[0], localMatrix[11]);
	globalMatrix.addElem(elem.info[6], elem.info[2], localMatrix[12]);
	globalMatrix.addElem(elem.info[6], elem.info[4], localMatrix[13]);
	globalMatrix.addElem(elem.info[6], elem.info[5], localMatrix[14]);
}

void Task::addLocalRigtPartToGlobal(const FiniteElem& elem) {
	/*
		!!! CHANGE IT !!!
	*/
	// TODO: consider new mapping local to global;
	f[elem.info[0]] += localRightPart[2];
	f[elem.info[1]] += localRightPart[1];
	f[elem.info[2]] += localRightPart[3];
	f[elem.info[3]] += localRightPart[0];
	f[elem.info[4]] += localRightPart[4];
	f[elem.info[5]] += localRightPart[5];
	f[elem.info[6]] += localRightPart[6];
}

void Task::setFirstBoundaryConditions() {
	/*
	!!! CHANGE IT !!!
*/

	// TODO: what we have with the Boundary conditions?;
	for (int ind : boundariesValue) {
		globalMatrix.setFirstBoundaryCondition(ind);
		f[ind] = 100;
		throw "not implemented";
	}
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


void Task::solve() {
	for (auto& el : elems) {
		calculateLocalMatrix(el);
		//calculateLocalRightPart(el);

		addLocalMatrixToGlobal(el);
		//addLocalRigtPartToGlobal(el);
	}

	// Дублирование элементов из нижнего треугольника матрицы в верхний
	globalMatrix.fillGGU();

	setFirstBoundaryConditions();
	PARDISOsolve();
}