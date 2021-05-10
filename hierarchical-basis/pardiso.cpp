#include "stdafx.h"
#include "pardiso.h"

extern ofstream logfile;

pardiso_solver::pardiso_solver()
{
	ia=NULL;
	ja=NULL;
	a=NULL;
	perm=NULL;

	n = 0;

	init();
}

pardiso_solver::~pardiso_solver()
{
	clear();
}

void pardiso_solver::init()
{
	int i;

	mtype = 11; //2
	nrhs = 1;
	maxfct = 1;
	mnum = 1;
	msglvl = 1;
	phase = 13;

	for(i=0;i<64;i++){pt[i]=0;}
	for(i=0;i<64;i++){iparm[i]=0;}

	msglvl=0;
}

void pardiso_solver::clear()
{
	if(a){delete [] a;  a=NULL;}
	if(ia){delete [] ia; ia=NULL;}
	if(ja){delete [] ja; ja=NULL;}
	if(perm){delete [] perm; perm=NULL;}
}

void pardiso_solver::factorize(int nb,int *ig,int *jg, double *ggl, double *ggu, double *di,int nthreads)
{
	int i;
	int ig_n_1=0;
	int sz_iptr=0;
	int sz_jptr=0;

	init();

	if(n!=nb)
	{
		clear();
	}

	mkl_set_num_threads(nthreads);

	// ��������������� � CSR
	f.From2x2ToCSR_Real_1_NS(nb, ig, &sz_iptr, &sz_jptr);

	if(n!=nb)
	{
		ia = new MKL_INT64[sz_iptr];
		ja = new MKL_INT64[sz_jptr];
		a = new double[sz_jptr];
	}

	f.FromRSFToCSR_Real_2_NS(nb, ig, jg, di, ggl, ggu, ia, ja, a);

	for(i=0;i<sz_iptr;i++){ia[i]++;}
	for(i=0;i<sz_jptr;i++){ja[i]++;}

	if(n!=nb)
	{
		perm = new MKL_INT64[nb];
	}

	phase = 12;
	n = nb;

	pardiso_64(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, perm, &nrhs, iparm, &msglvl, NULL, NULL, &info);
	cout << "pardiso factorize info = " << info << endl;
}

void pardiso_solver::solve_nrhs(int _nrhs,double *pr,double *x)
{
	int i;
	for(i=0;i<n;i++){x[i]=0.0;}
	phase = 33;
	nrhs = _nrhs;
	pardiso_64(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, perm, &nrhs, iparm, &msglvl, pr, x, &info);
	cout << "pardiso solve info = " << info << endl;
}

void pardiso_solver::stop_solver()
{
	phase = -1;
	pardiso_64(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, perm, &nrhs, iparm, &msglvl, NULL, NULL, &info);
}
