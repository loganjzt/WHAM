// 
//	getPvN.dat
// 
// Created by Zhitong on Nov 12 2018
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

int get_PvN()
{

	x1Max = max_int_2darray(nSimu,nDataMax,x1);
	x1Min = min_int_2darray(nSimu,nDataMax,x1);

	x1TildaMax = max_double_2darray(nSimu,nDataMax,x1Tilda);
	x1TildaMin = min_double_2darray(nSimu,nDataMax,x1Tilda);

	px1 = new double [x1Max - x1Min + 1];
	px1Tilda = new double [int( (x1TildaMax-x1TildaMin) / dx1Tilda + 1)];
	for(int i = 0 ; i < x1Max - x1Min + 1; i ++) px1[i] = 0.0;
	for(int i = 0 ; i < int( (x1TildaMax-x1TildaMin) / dx1Tilda + 1); i ++) px1Tilda[i] = 0.0;

	FILE* outputData;
	double normC;	// denominator

	double* sumC = new double [nSimu];
	double* sumC_d = new double [nDataMax*nSimu];

	int ii = 0;		// index for sumC_d
	// calculate normC , INDUS paper eq 7a
	for(int i = 0; i < nSimu; i++){
		for(int j = 0; j < nData[i] ; j++){
			for(int k = 0; k < nSimu ; k++){
				sumC[k] = log(nData[k])+logebWk[k][i][j]+logebfk[k];
			}		
			sumC_d[ii] = -logsumexp(nSimu,sumC);	
			ii += 1;
		}
    }

	// ii shoulde be equal to nData total 
    normC = -logsumexp(ii,sumC_d);
	printf("# normC is %g\n",normC);

	for(int i = 0; i < nSimu;i++){
		for(int j = 0; j < nData[i] ; j++){
			for(int k = 0; k < nSimu ; k++ ){
				sumC[k] = log(nData[k])+logebWk[k][i][j]+logebfk[k];
			}
			px1[int(x1[i][j] - x1Min )] += exp( normC - logsumexp(nSimu,sumC));
			px1Tilda[int((x1Tilda[i][j] - x1TildaMin ) / dx1Tilda ) ] += exp( normC - logsumexp(nSimu,sumC));
		}
    }

	//pvN
    outputData = fopen("pvN.dat","w");
    if(outputData == NULL){
		printf("# error in creating outputData file pvN.dat \n");
		return 1;
	}
	fprintf(outputData,"# n\tpvN\t-lnP\n");

    for(int i = 0; i <= (x1Max - x1Min) ; i++){
        fprintf(outputData,"%d\t%e\t%e\n", i + x1Min ,px1[i],-log(px1[i]));
    }

    fclose(outputData);

	//pvNTilda    
    outputData = fopen("pvNt.dat","w");
    if(outputData == NULL){
		printf("# error in creating outputData file pvNt.dat \n");
		return 1;
	}
	fprintf(outputData,"# n\tpvN\t-lnP\n");

    for(int i = 0; i <= int( (x1TildaMax - x1TildaMin)/dx1Tilda ); i++){
        fprintf(outputData,"%.2f\t%e\t%e\n",i + x1TildaMin, px1Tilda[i], -log(px1Tilda[i]));
    }

    fclose(outputData);


	delete [] sumC;
	delete [] sumC_d;

	return 0;
}
