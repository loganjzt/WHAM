//
// test data2ebWk.cpp
//
// Created by Zhitong Nov 12,2018
//
#include <stdio.h>
#include <cstdio>


	using namespace std;

int nSimu;
int *nData;
double ***logebWk;
double **x1Tilda;
double *x1Star;
double *mu, *kappa;

#include "/home/zhitongj/Git/wham/data2ebWk.cpp"

int main(){

	nSimu = 1;

	nData = new int [nSimu];
	nData[0] = 1;

	mu = new double [nSimu];
	kappa = new double [nSimu];
	x1Star = new double [nSimu];
	
	mu[0] = 1.0; 
	kappa[0] = 0.1;	
	x1Star[0] = 2.0;

	x1Tilda = new double* [nSimu];
	for(int i = 0; i < nSimu; i++){
		x1Tilda[i] = new double [nData[i]];
		for(int j = 0; j < nData[i] ; j++){
			x1Tilda[i][j] = 3.0;
		}
	}

	logebWk = new double** [nSimu];
	for(int i = 0; i < nSimu; i++){
		logebWk[i] = new double* [nSimu];
		for(int j = 0 ; j < nSimu ; j++){
			logebWk[i][j] = new double [nData[j]];
		}
	}

	double tmp = linear_umbrella(x1Tilda[0][0],mu[0],kappa[0],x1Star[0]);
	printf("# tmp = %g \n",tmp);

	data2ebWk(linear_umbrella);

	for(int i = 0; i < nSimu; i++){
		for(int j = 0 ; j < nSimu ; j++){
			for(int k = 0 ; k < nData[j]; k++){
				printf("# mu = %g , kappa = %g , nStar = %g , nTilda = %g , ebWk = %g \n",mu[i],kappa[i],x1Star[i],x1Tilda[j][k],logebWk[i][j][k]);
			}
		}
	}

	return 0;
}


