#pragma once
#include"stdafx.h"
#include"utils.h"

struct Matrix3D
{
	double M_R[3][3];
};

struct Table
{
	int sz;
	vector<double> arg, val;
	void Read(ifstream &inf)
	{
		int i;
		inf >> sz;
		arg.resize(sz);
		val.resize(sz);
		for (i = 0; i < sz; i++){ inf >> arg[i] >> val[i]; }
	}
	double GetValue(double x)
	{
		int i;
		double y, cff;
		y = 0.0;
		if (sz)
		{
			i = FindIntervalInDoubleVec(arg, x);
			if (i != -1)
			{
				cff = (x - arg[i]) / (arg[i + 1] - arg[i]);
				y = (1.0 - cff)*val[i] + cff*val[i + 1];
			}
			else if (x < arg[0]){ y = val[0]; }
			else { y = val[sz - 1]; }				/*(x>arg[sz-1])*/
		}
		return y;
	}
};

struct Point3D
{
	double crd[3];
};

struct Element3D
{
	static const int el_n = 8;
	int nodes[8], edges[12], faces[6], basenumer[64];
	int mat, type, sizeLocMatrix, type_curv, dim_curv;
};

struct Element3DForm
{
	//static int type_curv; // тип базиса, используемого для описания формы КЭ (2 - триквадратичные, 3 - трикубические)
	int el_n;
	int nodes[64];
};

struct Face
{
	int elem, gran;
};

struct SimpleSource : public Face
{
	double S;
};

struct Source : public Face
{
	double S[3];
};

struct KU_1
{
	int V[3];
	double U[3];
};

struct StiffProp
{
	static double A[6][6], S[6][6];

	bool fexist, fuseg;
	double Etcff, TEmin, TEmax;
	double E_TH[3], E[3], nu[3], G[3], h;
	double con[3], scon[3], ccon[3];
	double M_X[3][3], M_Y[3][3], M_Z[3][3], I_M_X[3][3], I_M_Y[3][3], I_M_Z[3][3], M_T[3][3], M_R[3][3], I_M_R[3][3];
	double E_R[6][6], S_R[6][6], T_LU[6][6], V_LU[6];

	double GetTEmin()
	{
		return TEmin;
	}

	void GetETH(double *pE_TH)
	{
		pE_TH[0] = E_TH[0];
		pE_TH[1] = E_TH[1];
		pE_TH[2] = E_TH[2];
	}

	void CalculateMatrixes()
	{
		int i;

		for (i = 0; i < 3; i++)
		{
			scon[i] = sin(con[i]);
			ccon[i] = cos(con[i]);
		}

		M_X[0][0] = 1.00000;		M_X[0][1] = 0.00000;		M_X[0][2] = 0.00000;
		M_X[1][0] = 0.00000;		M_X[1][1] = ccon[0];		M_X[1][2] = -scon[0];
		M_X[2][0] = 0.00000;		M_X[2][1] = scon[0];		M_X[2][2] = ccon[0];

		M_Y[0][0] = ccon[1];		M_Y[0][1] = 0.00000;		M_Y[0][2] = scon[1];
		M_Y[1][0] = 0.00000;		M_Y[1][1] = 1.00000;		M_Y[1][2] = 0.00000;
		M_Y[2][0] = -scon[1];		M_Y[2][1] = 0.00000;		M_Y[2][2] = ccon[1];

		M_Z[0][0] = ccon[2];		M_Z[0][1] = -scon[2];		M_Z[0][2] = 0.00000;
		M_Z[1][0] = scon[2];		M_Z[1][1] = ccon[2];		M_Z[1][2] = 0.00000;
		M_Z[2][0] = 0.00000;		M_Z[2][1] = 0.00000;		M_Z[2][2] = 1.00000;

		I_M_X[0][0] = 1.00000;	I_M_X[0][1] = 0.00000;	I_M_X[0][2] = 0.00000;
		I_M_X[1][0] = 0.00000;	I_M_X[1][1] = ccon[0];	I_M_X[1][2] = scon[0];
		I_M_X[2][0] = 0.00000;	I_M_X[2][1] = -scon[0];	I_M_X[2][2] = ccon[0];

		I_M_Y[0][0] = ccon[1];	I_M_Y[0][1] = 0.00000;	I_M_Y[0][2] = -scon[1];
		I_M_Y[1][0] = 0.00000;	I_M_Y[1][1] = 1.00000;	I_M_Y[1][2] = 0.00000;
		I_M_Y[2][0] = scon[1];	I_M_Y[2][1] = 0.00000;	I_M_Y[2][2] = ccon[1];

		I_M_Z[0][0] = ccon[2];	I_M_Z[0][1] = scon[2];	I_M_Z[0][2] = 0.00000;
		I_M_Z[1][0] = -scon[2];	I_M_Z[1][1] = ccon[2];	I_M_Z[1][2] = 0.00000;
		I_M_Z[2][0] = 0.00000;	I_M_Z[2][1] = 0.00000;	I_M_Z[2][2] = 1.00000;


		Mult_Matrix((double*)M_X, (double*)M_Y, (double*)M_T, 3);		// T=X*Y
		Mult_Matrix((double*)M_T, (double*)M_Z, (double*)M_R, 3);		// R=T*Z=X*Y*Z

		Mult_Matrix((double*)I_M_Z, (double*)I_M_Y, (double*)M_T, 3);	// T=I_Z*I_Y
		Mult_Matrix((double*)M_T, (double*)I_M_X, (double*)I_M_R, 3);	// I_R=T*I_X=I_Z*I_Y*I_X


		////Transpose((double *)M_R,3);
		////Transpose((double *)I_M_R,3);

		//// Ex,Ey,Ez
		//for(i=0;i<3;i++)
		//{
		//	// Ex,Ey,Ez
		//	E_R[i][0]=M_R[0][i]*M_R[0][i];
		//	E_R[i][1]=M_R[1][i]*M_R[1][i];
		//	E_R[i][2]=M_R[2][i]*M_R[2][i];
		//	// Exy,Eyz,Exz
		//	E_R[i][3]=M_R[0][i]*M_R[1][i];
		//	E_R[i][4]=M_R[1][i]*M_R[2][i];
		//	E_R[i][5]=M_R[0][i]*M_R[2][i];
		//}

		//// Exy: Ex,Ey,Ez
		//E_R[3][0]=2.0*M_R[0][0]*M_R[0][1];
		//E_R[3][1]=2.0*M_R[1][0]*M_R[1][1];
		//E_R[3][2]=2.0*M_R[2][0]*M_R[2][1];
		//// Exy: Exy,Eyz,Exz
		//E_R[3][3]=M_R[0][0]*M_R[1][1]+M_R[1][0]*M_R[0][1];
		//E_R[3][4]=M_R[1][0]*M_R[2][1]+M_R[2][0]*M_R[1][1];
		//E_R[3][5]=M_R[0][0]*M_R[2][1]+M_R[2][0]*M_R[0][1];

		//// Exz: Ex,Ey,Ez
		//E_R[4][0]=2.0*M_R[0][0]*M_R[0][2];
		//E_R[4][1]=2.0*M_R[1][0]*M_R[1][2];
		//E_R[4][2]=2.0*M_R[2][0]*M_R[2][2];
		//// Exz: Exy,Eyz,Exz
		//E_R[4][3]=M_R[0][0]*M_R[1][2]+M_R[1][0]*M_R[0][2];
		//E_R[4][4]=M_R[1][1]*M_R[2][2]+M_R[2][1]*M_R[1][2];
		//E_R[4][5]=M_R[0][1]*M_R[2][2]+M_R[2][1]*M_R[0][2];

		//// Eyz: Ex,Ey,Ez
		//E_R[5][0]=2.0*M_R[0][1]*M_R[0][2];
		//E_R[5][1]=2.0*M_R[1][1]*M_R[1][2];
		//E_R[5][2]=2.0*M_R[2][1]*M_R[2][2];
		//// Eyz: Exy,Eyz,Exz
		//E_R[5][3]=M_R[0][1]*M_R[1][2]+M_R[1][1]*M_R[0][2];
		//E_R[5][4]=M_R[1][0]*M_R[2][2]+M_R[2][0]*M_R[1][2];
		//E_R[5][5]=M_R[0][0]*M_R[2][2]+M_R[2][0]*M_R[0][2];

		//// Sx,Sy,Sz
		//for(i=0;i<3;i++)
		//{
		//	// Sx,Sy,Sz
		//	S_R[i][0]=M_R[0][i]*M_R[0][i];
		//	S_R[i][1]=M_R[1][i]*M_R[1][i];
		//	S_R[i][2]=M_R[2][i]*M_R[2][i];
		//	// Sxy,Syz,Sxz
		//	S_R[i][3]=2.0*M_R[0][i]*M_R[1][i];
		//	S_R[i][4]=2.0*M_R[1][i]*M_R[2][i];
		//	S_R[i][5]=2.0*M_R[0][i]*M_R[2][i];
		//}

		//// Sxy: Sx,Sy,Sz
		//S_R[3][0]=M_R[0][0]*M_R[0][1];
		//S_R[3][1]=M_R[1][0]*M_R[1][1];
		//S_R[3][2]=M_R[2][0]*M_R[2][1];
		//// Sxy: Sxy,Syz,Sxz
		//S_R[3][3]=M_R[0][0]*M_R[1][1]+M_R[1][0]*M_R[0][1];
		//S_R[3][4]=M_R[1][0]*M_R[2][1]+M_R[2][0]*M_R[1][1];
		//S_R[3][5]=M_R[0][0]*M_R[2][1]+M_R[2][0]*M_R[0][1];

		//// Sxz: Sx,Sy,Sz
		//S_R[4][0]=M_R[0][0]*M_R[0][2];
		//S_R[4][1]=M_R[1][0]*M_R[1][2];
		//S_R[4][2]=M_R[2][0]*M_R[2][2];
		//// Sxz: Sxy,Syz,Sxz
		//S_R[4][3]=M_R[0][0]*M_R[1][2]+M_R[1][0]*M_R[0][2];
		//S_R[4][4]=M_R[1][1]*M_R[2][2]+M_R[2][1]*M_R[1][2];
		//S_R[4][5]=M_R[0][1]*M_R[2][2]+M_R[2][1]*M_R[0][2];

		//// Syz: Sx,Sy,Sz
		//S_R[5][0]=M_R[0][1]*M_R[0][2];
		//S_R[5][1]=M_R[1][1]*M_R[1][2];
		//S_R[5][2]=M_R[2][1]*M_R[2][2];
		//// Syz: Sxy,Syz,Sxz
		//S_R[5][3]=M_R[0][1]*M_R[1][2]+M_R[1][1]*M_R[0][2];
		//S_R[5][4]=M_R[1][0]*M_R[2][2]+M_R[2][0]*M_R[1][2];
		//S_R[5][5]=M_R[0][0]*M_R[2][2]+M_R[2][0]*M_R[0][2];


		// Ex,Ey,Ez
		for (i = 0; i < 3; i++)
		{
			// Ex,Ey,Ez
			E_R[i][0] = M_R[i][0] * M_R[i][0];
			E_R[i][1] = M_R[i][1] * M_R[i][1];
			E_R[i][2] = M_R[i][2] * M_R[i][2];
			// Exy,Eyz,Exz
			E_R[i][3] = M_R[i][0] * M_R[i][1];
			E_R[i][4] = M_R[i][1] * M_R[i][2];
			E_R[i][5] = M_R[i][0] * M_R[i][2];
		}

		// Exy: Ex,Ey,Ez
		E_R[3][0] = 2.0*M_R[0][0] * M_R[1][0];
		E_R[3][1] = 2.0*M_R[0][1] * M_R[1][1];
		E_R[3][2] = 2.0*M_R[0][2] * M_R[1][2];
		// Exy: Exy,Eyz,Exz
		E_R[3][3] = M_R[0][0] * M_R[1][1] + M_R[0][1] * M_R[1][0];
		E_R[3][4] = M_R[0][1] * M_R[1][2] + M_R[0][2] * M_R[1][1];
		E_R[3][5] = M_R[0][0] * M_R[1][2] + M_R[0][2] * M_R[1][0];

		// Exz: Ex,Ey,Ez
		E_R[4][0] = 2.0*M_R[1][0] * M_R[2][0];
		E_R[4][1] = 2.0*M_R[1][1] * M_R[2][1];
		E_R[4][2] = 2.0*M_R[1][2] * M_R[2][2];
		// Exz: Exy,Eyz,Exz
		E_R[4][3] = M_R[1][0] * M_R[2][1] + M_R[1][1] * M_R[2][0];
		E_R[4][4] = M_R[1][1] * M_R[2][2] + M_R[1][2] * M_R[2][1];
		E_R[4][5] = M_R[1][0] * M_R[2][2] + M_R[1][2] * M_R[2][0];

		// Eyz: Ex,Ey,Ez
		E_R[5][0] = 2.0*M_R[0][0] * M_R[2][0];
		E_R[5][1] = 2.0*M_R[0][1] * M_R[2][1];
		E_R[5][2] = 2.0*M_R[0][2] * M_R[2][2];
		// Eyz: Exy,Eyz,Exz
		E_R[5][3] = M_R[0][0] * M_R[2][1] + M_R[0][1] * M_R[2][0];
		E_R[5][4] = M_R[0][1] * M_R[2][2] + M_R[0][2] * M_R[2][1];
		E_R[5][5] = M_R[0][0] * M_R[2][2] + M_R[0][2] * M_R[2][0];

		// Sx,Sy,Sz
		for (i = 0; i < 3; i++)
		{
			// Sx,Sy,Sz
			S_R[i][0] = M_R[i][0] * M_R[i][0];
			S_R[i][1] = M_R[i][1] * M_R[i][1];
			S_R[i][2] = M_R[i][2] * M_R[i][2];
			// Sxy,Syz,Sxz
			S_R[i][3] = 2.0*M_R[i][0] * M_R[i][1];
			S_R[i][4] = 2.0*M_R[i][1] * M_R[i][2];
			S_R[i][5] = 2.0*M_R[i][0] * M_R[i][2];
		}

		// Sxy: Sx,Sy,Sz
		S_R[3][0] = M_R[0][0] * M_R[1][0];
		S_R[3][1] = M_R[0][1] * M_R[1][1];
		S_R[3][2] = M_R[0][2] * M_R[1][2];
		// Sxy: Sxy,Syz,Sxz
		S_R[3][3] = M_R[0][0] * M_R[1][1] + M_R[0][1] * M_R[1][0];
		S_R[3][4] = M_R[0][1] * M_R[1][2] + M_R[0][2] * M_R[1][1];
		S_R[3][5] = M_R[0][0] * M_R[1][2] + M_R[0][2] * M_R[1][0];

		// Sxz: Sx,Sy,Sz
		S_R[4][0] = M_R[1][0] * M_R[2][0];
		S_R[4][1] = M_R[1][1] * M_R[2][1];
		S_R[4][2] = M_R[1][2] * M_R[2][2];
		// Sxz: Sxy,Syz,Sxz
		S_R[4][3] = M_R[1][0] * M_R[2][1] + M_R[1][1] * M_R[2][0];
		S_R[4][4] = M_R[1][1] * M_R[2][2] + M_R[1][2] * M_R[2][1];
		S_R[4][5] = M_R[1][0] * M_R[2][2] + M_R[1][2] * M_R[2][0];

		// Syz: Sx,Sy,Sz
		S_R[5][0] = M_R[0][0] * M_R[2][0];
		S_R[5][1] = M_R[0][1] * M_R[2][1];
		S_R[5][2] = M_R[0][2] * M_R[2][2];
		// Syz: Sxy,Syz,Sxz
		S_R[5][3] = M_R[0][0] * M_R[2][1] + M_R[0][1] * M_R[2][0];
		S_R[5][4] = M_R[0][1] * M_R[2][2] + M_R[0][2] * M_R[2][1];
		S_R[5][5] = M_R[0][0] * M_R[2][2] + M_R[0][2] * M_R[2][0];
	}

	void CalculateMatrixes(double p_M_R[][3])
	{
		int i, j;

		for (i = 0; i < 3; i++)
		{
			for (j = 0; j < 3; j++)
			{
				I_M_R[i][j] = M_R[i][j] = p_M_R[i][j];
			}
		}
		Transpose((double *)I_M_R, 3);

		////Transpose((double *)M_R,3);
		////Transpose((double *)I_M_R,3);

		//// Ex,Ey,Ez
		//for(i=0;i<3;i++)
		//{
		//	// Ex,Ey,Ez
		//	E_R[i][0]=M_R[0][i]*M_R[0][i];
		//	E_R[i][1]=M_R[1][i]*M_R[1][i];
		//	E_R[i][2]=M_R[2][i]*M_R[2][i];
		//	// Exy,Eyz,Exz
		//	E_R[i][3]=M_R[0][i]*M_R[1][i];
		//	E_R[i][4]=M_R[1][i]*M_R[2][i];
		//	E_R[i][5]=M_R[0][i]*M_R[2][i];
		//}

		//// Exy: Ex,Ey,Ez
		//E_R[3][0]=2.0*M_R[0][0]*M_R[0][1];
		//E_R[3][1]=2.0*M_R[1][0]*M_R[1][1];
		//E_R[3][2]=2.0*M_R[2][0]*M_R[2][1];
		//// Exy: Exy,Eyz,Exz
		//E_R[3][3]=M_R[0][0]*M_R[1][1]+M_R[1][0]*M_R[0][1];
		//E_R[3][4]=M_R[1][0]*M_R[2][1]+M_R[2][0]*M_R[1][1];
		//E_R[3][5]=M_R[0][0]*M_R[2][1]+M_R[2][0]*M_R[0][1];

		//// Exz: Ex,Ey,Ez
		//E_R[4][0]=2.0*M_R[0][0]*M_R[0][2];
		//E_R[4][1]=2.0*M_R[1][0]*M_R[1][2];
		//E_R[4][2]=2.0*M_R[2][0]*M_R[2][2];
		//// Exz: Exy,Eyz,Exz
		//E_R[4][3]=M_R[0][0]*M_R[1][2]+M_R[1][0]*M_R[0][2];
		//E_R[4][4]=M_R[1][1]*M_R[2][2]+M_R[2][1]*M_R[1][2];
		//E_R[4][5]=M_R[0][1]*M_R[2][2]+M_R[2][1]*M_R[0][2];

		//// Eyz: Ex,Ey,Ez
		//E_R[5][0]=2.0*M_R[0][1]*M_R[0][2];
		//E_R[5][1]=2.0*M_R[1][1]*M_R[1][2];
		//E_R[5][2]=2.0*M_R[2][1]*M_R[2][2];
		//// Eyz: Exy,Eyz,Exz
		//E_R[5][3]=M_R[0][1]*M_R[1][2]+M_R[1][1]*M_R[0][2];
		//E_R[5][4]=M_R[1][0]*M_R[2][2]+M_R[2][0]*M_R[1][2];
		//E_R[5][5]=M_R[0][0]*M_R[2][2]+M_R[2][0]*M_R[0][2];

		//// Sx,Sy,Sz
		//for(i=0;i<3;i++)
		//{
		//	// Sx,Sy,Sz
		//	S_R[i][0]=M_R[0][i]*M_R[0][i];
		//	S_R[i][1]=M_R[1][i]*M_R[1][i];
		//	S_R[i][2]=M_R[2][i]*M_R[2][i];
		//	// Sxy,Syz,Sxz
		//	S_R[i][3]=2.0*M_R[0][i]*M_R[1][i];
		//	S_R[i][4]=2.0*M_R[1][i]*M_R[2][i];
		//	S_R[i][5]=2.0*M_R[0][i]*M_R[2][i];
		//}

		//// Sxy: Sx,Sy,Sz
		//S_R[3][0]=M_R[0][0]*M_R[0][1];
		//S_R[3][1]=M_R[1][0]*M_R[1][1];
		//S_R[3][2]=M_R[2][0]*M_R[2][1];
		//// Sxy: Sxy,Syz,Sxz
		//S_R[3][3]=M_R[0][0]*M_R[1][1]+M_R[1][0]*M_R[0][1];
		//S_R[3][4]=M_R[1][0]*M_R[2][1]+M_R[2][0]*M_R[1][1];
		//S_R[3][5]=M_R[0][0]*M_R[2][1]+M_R[2][0]*M_R[0][1];

		//// Sxz: Sx,Sy,Sz
		//S_R[4][0]=M_R[0][0]*M_R[0][2];
		//S_R[4][1]=M_R[1][0]*M_R[1][2];
		//S_R[4][2]=M_R[2][0]*M_R[2][2];
		//// Sxz: Sxy,Syz,Sxz
		//S_R[4][3]=M_R[0][0]*M_R[1][2]+M_R[1][0]*M_R[0][2];
		//S_R[4][4]=M_R[1][1]*M_R[2][2]+M_R[2][1]*M_R[1][2];
		//S_R[4][5]=M_R[0][1]*M_R[2][2]+M_R[2][1]*M_R[0][2];

		//// Syz: Sx,Sy,Sz
		//S_R[5][0]=M_R[0][1]*M_R[0][2];
		//S_R[5][1]=M_R[1][1]*M_R[1][2];
		//S_R[5][2]=M_R[2][1]*M_R[2][2];
		//// Syz: Sxy,Syz,Sxz
		//S_R[5][3]=M_R[0][1]*M_R[1][2]+M_R[1][1]*M_R[0][2];
		//S_R[5][4]=M_R[1][0]*M_R[2][2]+M_R[2][0]*M_R[1][2];
		//S_R[5][5]=M_R[0][0]*M_R[2][2]+M_R[2][0]*M_R[0][2];


		// Ex,Ey,Ez
		for (i = 0; i < 3; i++)
		{
			// Ex,Ey,Ez
			E_R[i][0] = M_R[i][0] * M_R[i][0];
			E_R[i][1] = M_R[i][1] * M_R[i][1];
			E_R[i][2] = M_R[i][2] * M_R[i][2];
			// Exy,Eyz,Exz
			E_R[i][3] = M_R[i][0] * M_R[i][1];
			E_R[i][4] = M_R[i][1] * M_R[i][2];
			E_R[i][5] = M_R[i][0] * M_R[i][2];
		}

		// Exy: Ex,Ey,Ez
		E_R[3][0] = 2.0*M_R[0][0] * M_R[1][0];
		E_R[3][1] = 2.0*M_R[0][1] * M_R[1][1];
		E_R[3][2] = 2.0*M_R[0][2] * M_R[1][2];
		// Exy: Exy,Eyz,Exz
		E_R[3][3] = M_R[0][0] * M_R[1][1] + M_R[0][1] * M_R[1][0];
		E_R[3][4] = M_R[0][1] * M_R[1][2] + M_R[0][2] * M_R[1][1];
		E_R[3][5] = M_R[0][0] * M_R[1][2] + M_R[0][2] * M_R[1][0];

		// Exz: Ex,Ey,Ez
		E_R[4][0] = 2.0*M_R[1][0] * M_R[2][0];
		E_R[4][1] = 2.0*M_R[1][1] * M_R[2][1];
		E_R[4][2] = 2.0*M_R[1][2] * M_R[2][2];
		// Exz: Exy,Eyz,Exz
		E_R[4][3] = M_R[1][0] * M_R[2][1] + M_R[1][1] * M_R[2][0];
		E_R[4][4] = M_R[1][1] * M_R[2][2] + M_R[1][2] * M_R[2][1];
		E_R[4][5] = M_R[1][0] * M_R[2][2] + M_R[1][2] * M_R[2][0];

		// Eyz: Ex,Ey,Ez
		E_R[5][0] = 2.0*M_R[0][0] * M_R[2][0];
		E_R[5][1] = 2.0*M_R[0][1] * M_R[2][1];
		E_R[5][2] = 2.0*M_R[0][2] * M_R[2][2];
		// Eyz: Exy,Eyz,Exz
		E_R[5][3] = M_R[0][0] * M_R[2][1] + M_R[0][1] * M_R[2][0];
		E_R[5][4] = M_R[0][1] * M_R[2][2] + M_R[0][2] * M_R[2][1];
		E_R[5][5] = M_R[0][0] * M_R[2][2] + M_R[0][2] * M_R[2][0];

		// Sx,Sy,Sz
		for (i = 0; i < 3; i++)
		{
			// Sx,Sy,Sz
			S_R[i][0] = M_R[i][0] * M_R[i][0];
			S_R[i][1] = M_R[i][1] * M_R[i][1];
			S_R[i][2] = M_R[i][2] * M_R[i][2];
			// Sxy,Syz,Sxz
			S_R[i][3] = 2.0*M_R[i][0] * M_R[i][1];
			S_R[i][4] = 2.0*M_R[i][1] * M_R[i][2];
			S_R[i][5] = 2.0*M_R[i][0] * M_R[i][2];
		}

		// Sxy: Sx,Sy,Sz
		S_R[3][0] = M_R[0][0] * M_R[1][0];
		S_R[3][1] = M_R[0][1] * M_R[1][1];
		S_R[3][2] = M_R[0][2] * M_R[1][2];
		// Sxy: Sxy,Syz,Sxz
		S_R[3][3] = M_R[0][0] * M_R[1][1] + M_R[0][1] * M_R[1][0];
		S_R[3][4] = M_R[0][1] * M_R[1][2] + M_R[0][2] * M_R[1][1];
		S_R[3][5] = M_R[0][0] * M_R[1][2] + M_R[0][2] * M_R[1][0];

		// Sxz: Sx,Sy,Sz
		S_R[4][0] = M_R[1][0] * M_R[2][0];
		S_R[4][1] = M_R[1][1] * M_R[2][1];
		S_R[4][2] = M_R[1][2] * M_R[2][2];
		// Sxz: Sxy,Syz,Sxz
		S_R[4][3] = M_R[1][0] * M_R[2][1] + M_R[1][1] * M_R[2][0];
		S_R[4][4] = M_R[1][1] * M_R[2][2] + M_R[1][2] * M_R[2][1];
		S_R[4][5] = M_R[1][0] * M_R[2][2] + M_R[1][2] * M_R[2][0];

		// Syz: Sx,Sy,Sz
		S_R[5][0] = M_R[0][0] * M_R[2][0];
		S_R[5][1] = M_R[0][1] * M_R[2][1];
		S_R[5][2] = M_R[0][2] * M_R[2][2];
		// Syz: Sxy,Syz,Sxz
		S_R[5][3] = M_R[0][0] * M_R[2][1] + M_R[0][1] * M_R[2][0];
		S_R[5][4] = M_R[0][1] * M_R[2][2] + M_R[0][2] * M_R[2][1];
		S_R[5][5] = M_R[0][0] * M_R[2][2] + M_R[0][2] * M_R[2][0];
	}

	void CalculateMatrixD(double dT, double Tave, double(*D)[6], vector<Table> &YunTable, vector<Table> &PsnTable, bool f_calc_d = true)
	{
		int i, j;
		double Es[3], nus[3];

		//Transpose((double *)S_R,6);
		//Mult_Matrix((double*)S_R,(double*)E_R,(double*)D,6);

		for (i = 0; i < 3; i++)
		{
			Es[i] = E[i];
			nus[i] = nu[i];
		}

		//E[0]*=(1.0+Etcff*dT/(TEmax-TEmin));
		//E[1]*=(1.0+Etcff*dT/(TEmax-TEmin));
		//E[2]*=(1.0+Etcff*dT/(TEmax-TEmin));

		for (i = 0; i < 3; i++)
		{
			if (E[i] > 0.0 && YunTable[i].sz){ E[i] = YunTable[i].GetValue(Tave); }
			if (nu[i] > 0.0 && PsnTable[i].sz){ nu[i] = PsnTable[i].GetValue(Tave); }
		}

		if (f_calc_d)
		{

			h = 1.0 - nu[0] * nu[0] * (E[1] / E[0]) - nu[2] * nu[2] * (E[2] / E[1]) - nu[1] * nu[1] * (E[2] / E[0]) - 2 * nu[0] * nu[1] * nu[2] * (E[2] / E[0]);

			for (i = 0; i < 6; i++){ for (j = 0; j < 6; j++){ D[i][j] = 0.0; } }

			D[0][0] = (E[0] / h)*(1.0 - nu[2] * nu[2] * (E[2] / E[1]));
			D[1][1] = (E[1] / h)*(1.0 - nu[1] * nu[1] * (E[2] / E[0]));
			D[2][2] = (E[2] / h)*(1.0 - nu[0] * nu[0] * (E[1] / E[0]));

			D[0][1] = D[1][0] = (E[1] / h)*(nu[0] + nu[1] * nu[2] * (E[2] / E[1]));
			D[0][2] = D[2][0] = (E[2] / h)*(nu[1] + nu[2] * nu[0]);
			D[1][2] = D[2][1] = (E[2] / h)*(nu[2] + nu[1] * nu[0] * (E[1] / E[0]));

			if (fuseg)
			{
				D[3][3] = G[0];
				D[4][4] = G[2];
				D[5][5] = G[1];
			}
			else
			{
				D[3][3] = E[0] * E[1] / (E[0] * (1.0 + 2.0*nu[0] * E[1] / E[0]) + E[1]);
				D[4][4] = E[1] * E[2] / (E[1] * (1.0 + 2.0*nu[2] * E[2] / E[1]) + E[2]);
				D[5][5] = E[0] * E[2] / (E[0] * (1.0 + 2.0*nu[1] * E[2] / E[0]) + E[2]);
			}

		}

		Mult_Inverse_Gauss_Modify_Full((double *)S_R, (double *)D, (double *)T_LU, 6, (double *)A, (double *)S);
		Mult_Matrix((double*)T_LU, (double*)E_R, (double*)D, 6);

		for (i = 0; i < 3; i++)
		{
			E[i] = Es[i];
			nu[i] = nus[i];
		}

		//ofstream ofp;
		//ofp.open("D_mat",ios::app);
		//ofp<<scientific<<setprecision(14);
		//for(i=0;i<6;i++)
		//{
		//	for(j=0;j<6;j++)
		//	{
		//		ofp<<D[i][j]<<' ';
		//	}
		//	ofp<<'\n';
		//}
		//ofp<<'\n'<<'\n';
		//ofp.close();
		//ofp.clear();
	}
};
