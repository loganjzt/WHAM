//
// test getPvN.cpp
//
// Created by Zhitong Nov 12,2018
//
#include <stdio.h>
#include <cstdio>


	using namespace std;

int nSimu;
int *nData;
int nDataMax;
double ***logebWk;

double *logebfk;
int **x1;
int x1Max,x1Min;
double **x1Tilda;
double x1TildaMax,x1TildaMin;
double dx1Tilda;

double *px1,*px1Tilda;

#include "/home/zhitongj/Git/wham/getPvN.cpp"

int main(){

	nSimu = 2;

	dx1Tilda = 1.0;

	nData = new int [nSimu];

	for(int i = 0 ; i< nSimu; i++){
		nData[i] = 2;
	}	
	nDataMax = 2;

	x1 = new int* [nSimu];
	x1Tilda = new double* [nSimu];
	for(int i = 0; i < nSimu; i++){
		x1[i] = new int [nData[i]];
		x1Tilda[i] = new double [nData[i]];
	}

	logebfk = new double [nSimu];
	for(int i = 0; i < nSimu; i++) logebfk[i] = 1.0;
	
	logebWk = new double** [nSimu];
	for(int i = 0; i < nSimu ; i ++){
		logebWk[i] = new double* [nSimu]; 
		for(int j = 0; j < nSimu; j++){
			logebWk[i][j] = new double [nData[j]];
			for(int k = 0; k < nData[j]; k++){
				logebWk[i][j][k] = 0.0;
			}
		}
	}	


	x1[0][0] = 0.0;
	x1[0][1] = 1.0;
	x1[1][0] = 2.0;
	x1[1][1] = 1.0;

	x1Tilda[0][0] = 0.5;
	x1Tilda[0][1] = 1.5;
	x1Tilda[1][0] = 2.5;
	x1Tilda[1][1] = 1.5;

	get_PvN();


	return 0;
}


