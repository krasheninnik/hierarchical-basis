#pragma once
#include "FormatConverter.h"

struct pardiso_solver
{
	FormatConverter f;

	// матрица, сконвертированная в разреженный строчный формат
	MKL_INT64 *ia;
	MKL_INT64 *ja;
	double *a;

	// для PARDISO
	MKL_INT64 n;
	MKL_INT64 mtype;
	MKL_INT64 nrhs;
	void *pt[64];
	MKL_INT64 maxfct;
	MKL_INT64 mnum;
	MKL_INT64 msglvl;
	MKL_INT64 phase;
	MKL_INT64 *perm;
	MKL_INT64 iparm[64];
	MKL_INT64 info;

	pardiso_solver();
	~pardiso_solver();

	void factorize(int nb, int* ig, int* jg, double* ggl, double* ggu, double* di, int nthreads);
	void solve_nrhs(int nrhs,double *pr,double *x);
	void stop_solver();
	void clear();
	void init();
};
