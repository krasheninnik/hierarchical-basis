#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <direct.h>
#include <math.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <array>
#include <vector>
#include <set>
#include <complex>
#include <algorithm>

#define PI 3.1415926535897932

#define RETCODE_OK 0
#define RETCODE_NOMEM 1
#define RETCODE_NOFILE 2
#define RETCODE_OUTOFRANGE 3
#define RETCODE_ERROR 10

const double d_eps=1e-9;

const int size_i=sizeof(int);
const int size_d=sizeof(double);

using namespace std;

__forceinline ostream& operator < (ostream& file,const double& data)
{
	file.write((char*)&data,sizeof(data));
	return file;
}
__forceinline ostream& operator < (ostream& file,const int& data)
{
	file.write((char*)&data,sizeof(data));
	return file;
}

__forceinline istream& operator > (istream& file,double&  data)
{
	file.read((char*)&data,sizeof(data));
	return file;
}

__forceinline istream& operator > (istream& file,int&  data)
{
	file.read((char*)&data,sizeof(data));
	return file;
}

#include"mkl.h"
