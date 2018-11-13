//
//	test for save_whaminfo.cpp
//
#include <fstream>
#include <stdio.h>
#include <iostream>
#include <string>
#include <cstring>
#include <stdlib.h>
#include <cstdio> // fopen
#include "/home/zhitongj/Git/wham/misc_fn.cpp"

	using namespace std;

// test save_whaminfo fn
// variables declaration
int nSimu;
int *nSkip, *nData;
double *mu, *kappa, *x1Star, *logebfk;
char **fnarray;
double kBT;
double temperature;

// include code
#include "/home/zhitongj/Git/wham/save_whaminfo.cpp"

int main(){

	nSimu = 2;
	nSkip = new int [nSimu];
	nData = new int [nSimu];
	mu = new double [nSimu];
	kappa = new double [nSimu];	
	x1Star = new double [nSimu];
	logebfk = new double [nSimu];
	fnarray = Create2dCharArray(nSimu,500);

	kBT = 1.0;
	for(int i = 0 ; i < nSimu ; i++){
		nSkip[i] = 100;
		nData[i] = 200;
		mu[i] = 0.1;
		kappa[i] = 1.0;
		x1Star[i] = 33.0;
		logebfk[i] = 0.0;
		strcpy(fnarray[i] , "test.dat");
	}
	temperature = 300.0;
	
	printf("# Finish define variables\n");
	char* infoFileName = new char [500];
	strcpy(infoFileName,"test.whaminfo");

	save_whaminfo(infoFileName);

	return 0;
}
  

