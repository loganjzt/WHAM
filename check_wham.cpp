//
// check_wham.cpp
//
// calculate shan entropy to check wham
//
#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <cmath>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "/home/zhitongj/Git/wham/misc_fn.cpp" // logsumexp and min max fn

	using namespace std;

int check_wham(
	double (*potentialfn)(double, double, double, double), // potential fn
	int isave	// save reweighted distribution 
){


	double* s_shan = new double [nSimu];	// shan entropy of each window
	
	double** px1Tilda_sim = new double* [nSimu];	// distribution in each biased window
	double** logpx1Tilda_wham = new double* [nSimu];	// distribtuion from reweighting wham results
	double* lognormC_wham = new double [nSimu];

	int nx1TildaBin = int( ( x1TildaMax - x1TildaMin ) / dx1Tilda + 0.5 ) ; // number
	printf("# nx1TildaBin = %d ....\n",nx1TildaBin);
	for(int i = 0; i < nSimu; i++){
		px1Tilda_sim[i] = new double [nx1TildaBin];
		logpx1Tilda_wham[i] = new double [nx1TildaBin];	
	}


	for(int i = 0; i < nSimu; i++){
		s_shan[i] = 0.0;

		// bin simulation data
		for(int j = 0 ; j < nData[i] ; j++)
			px1Tilda_sim[i][int( (x1Tilda[i][j] - x1TildaMin )/dx1Tilda ) ] += 1.0;

		for(int j = 0; j < nx1TildaBin; j++  ){
			px1Tilda_sim[i][j] /= nData[i];		// normalize histograms from simulation
			logpx1Tilda_wham[i][j] = log(px1Tilda[j]) - (*potentialfn)(j*dx1Tilda+x1TildaMin,mu[i],kappa[i],x1Star[i]);
		}
		// nomalize histrogram from reweighting wham results
		lognormC_wham[i] = - 1.0 * logsumexp(nx1TildaBin, logpx1Tilda_wham[i]);
		for(int j = 0; j < nx1TildaBin; j++) logpx1Tilda_wham[i][j] += lognormC_wham[i];

		// calculate shan entropy
		for(int j = 0; j < nx1TildaBin; j++){
			if( px1Tilda_sim[i][j] > 1.0 / double(nData[i]) ){
				s_shan[i] += px1Tilda_sim[i][j] * abs( log(px1Tilda_sim[i][j]) - logpx1Tilda_wham[i][j] ) ;
			}
		}
	}
	
	if(isave == 1){
		printf("# saving entropy data...\n");
		FILE* checkFile = fopen("shanE.dat","w");
	
		fprintf(checkFile,"# kappa(kJ/mol)  Nstar    mu(kJ/mol)  S_shan\n");
		for(int i = 0; i < nSimu; i++){
			fprintf(checkFile,"%d\t%f\t%f\t%f\t%e\n",i+1,kappa[i],x1Star[i],mu[i],s_shan[i]);
		}
	}

	return 0;
}
