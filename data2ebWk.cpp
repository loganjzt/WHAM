//
// data2ebWk.cpp
//
// Created by Zhitong Nov 12,2018
//
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <cstdio>

	using namespace std;

int data2ebWk(double(*potentialfn)(double,double,double,double))
{

	logebWk = new double** [nSimu];
	for(int i = 0; i < nSimu; i++){
		logebWk[i] = new double* [nSimu];
		for(int j = 0 ; j < nSimu ; j++){
			logebWk[i][j] = new double [nDataMax];
		}
	}


	// calculate ebWk
	for(int i = 0; i < nSimu; i++){
		for(int j = 0; j < nData[i]; j++){
			for(int k = 0; k < nSimu ; k ++){
				logebWk[k][i][j] = (- 1.0 * (*potentialfn)(x1Tilda[i][j],mu[k],kappa[k],x1Star[k]) ) ;
			}
		}
	}

	return 0;
}


// Potentials
double linear_umbrella(double n, double m, double k, double nc)
{
   return(  m * n  +  0.5 * k * (n - nc) * (n - nc)  );
}

