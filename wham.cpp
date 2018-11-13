// 1D wham first version for Git, Apr 24 2017
//
// Modified on 11/09/2018

#define LBFGS_FLOAT 64
#include <stdio.h>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <cstring>
#include <string>
#include <math.h>
#include <fstream>
#include <cstdio>

	using namespace std;

// define variables
double temperature, kBT;

int nSimu;			// number of windows
char **dataFile;	// name of the input data files
int *nSkip;			// number of skiped data for each window
int *nData;			// number of loaded data for each window
double *lognData;	// log(nData[i])
int nDataMax;

double *mu,*x1Star,*kappa;	// parameter of biasing potential for each window

char **fnarray;
int **x1;
double **x1Tilda;	// x[nSimu][nData] and xTilda[nSimu][nData] time series 
int x1Max,x1Min;
double x1TildaMin,x1TildaMax;	// range of x and xTilda
double dx1,dx1Tilda;			// bin size of x1 and x1Tilda for histogram

double ***logebWk;		// [nSimu][nSimu][nData]	, ebWK(k,i,j) will store exp(-beta*U_k(nTilda(i,j)))
double *logebfk;		// free energy from adding biasing potential, used as intial values for lbfgs
double *difflogebfk;

double *px1,*px1Tilda,*logpx1Tilda;

#include "/home/zhitongj/Git/wham/loadData.cpp"
#include "/home/zhitongj/Git/wham/data2ebWk.cpp"
#include "/home/zhitongj/Git/wham/call_bfgs.cpp"
#include "/home/zhitongj/Git/wham/save_whaminfo.cpp"

#include "/home/zhitongj/Git/wham/getPvN.cpp"
#include "/home/zhitongj/Git/wham/check_wham.cpp"

int main(int argc, char *argv[])
{

	dx1Tilda = 1.0;

	if(argc < 1){
		printf("# wham info file name\n");
		return 1;
	}

	loadData(argv[1]);

	data2ebWk(linear_umbrella);
	call_bfgs();

	save_whaminfo(argv[1]);

	printf("# getting distribution...\n");
	get_PvN();

	printf("# checking shan entropy ...\n");
	check_wham(linear_umbrella,1);

    return 0;
}


