// WHAM code checked on Oct11 2016

#define LBFGS_FLOAT 64
#include <stdio.h>
#include "/home/pamish/packages/utils/lbfgs_install/include/lbfgs.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <cstring>
#include <string>
#include <math.h>
#include <fstream>
#include <malloc.h>
#include <cstdio>

	using namespace std;

    static int nFrame = 40000;
	static int nCut   = 10000;

double logsumexp(int nx, void *y);

double evaluate(void *instance,const double *x,double *g,const int n,const double step);

static int progress(void *instance,const double *x,const double *g,const double fx,
	const double xnorm,const double gnorm,const double step,int n,int k,int ls);

int main(int argc, char *argv[])
{
	int i,j,k = 0;

	const int nSimu = 19;				// Number of total Simulations
	const int nMin = -20;			// minimum Nstar
	const int dN = 10;

	/*BFGS input variables*/
	int ret = 0;
	double fx;
	double *x;
    lbfgs_parameter_t param;
 
    /*used to calculate Unmbrella potential */
	double temperature = 300.0;
    double beta = 1.0 / (temperature * 8.314 / 1000.0);

    /* get data from files and cal U */
	char f[40];
	FILE * data;
	
	double *kappa;
	int *nStar;
	double **nTw;
	int **ni;
	double time,n,ntw;
	int nMax = 0,nTmax = 0;		

	double ***umb;					//exp(-beta * umbrella potential for each window w nTw[j][k] 
	double *ci,C;					//free energy caused by adding biasing potential to each window;
	double *sumC_d, *sumC;			//sumC is for calculate C
	double *pvNt, *pvN;
    int *nTotal;

    /*allocate memory*/
	sumC = (double *)malloc(sizeof(double)*nSimu);
	sumC_d = (double *)malloc(sizeof(double)*(nSimu*nFrame));

	x = (double *)malloc(sizeof(double)*(nSimu-1));		// x[N-1]
	ci = (double *)malloc(sizeof(double)*nSimu);

	nStar = (int *)malloc(sizeof(int)*nSimu);
    kappa = (double *)malloc(sizeof(double)*nSimu);
	umb = (double***)malloc(sizeof(double **)*nSimu);

	ni = (int **)malloc(sizeof(int *)*nSimu);
	nTw = (double**)malloc(sizeof(double *)*nSimu);
	for(i = 0;i<nSimu;i++){
		ni[i] = (int *)malloc(sizeof(int)*nFrame);
		nTw[i] = (double *)malloc(sizeof(double)*nFrame);
		umb[i] = (double **)malloc(sizeof(double **)*nSimu);
		for(j = 0;j<nSimu;j++){umb[i][j] = (double *)malloc(sizeof(double)*nFrame);}
	}
	
	/* set nStar and kappa */
	for(i = 0;i<nSimu ;i++){kappa[i] = 0.500;}
	for(i = 0;i<nSimu;i++){nStar[i] = i*dN + nMin;}

	/*load files*/
	for(i = 0;i<nSimu ;i++){
		sprintf(f,"n%d/dynN.dat",nStar[i]);					// data FILES names
		//sprintf(f,"K%.3fN%d/dynN.dat",kappa[i],nStar[i]);
		data = fopen(f,"r");
		//cout<<f<<endl;
		if (data == NULL) printf("# Error in loading biased data, Nstar = %d\n",nStar[i]);
	
		for(j=0;j<nCut+5;j++){
			fgets(f , 100 , data);
		}

		//k = 0;
		//while(!feof(data)){
		
		for(k=0;k<nFrame;k++){
			fscanf(data,"%lf %lf %lf\n",&time,&n,&ntw);
			ni[i][k] = int(n+0.5);
			nTw[i][k] = ntw;
			//k++;
	    	if(int(ntw+0.5)>nTmax) nTmax = int(ntw+0.5);
			if(int(n+0.5)>nMax) nMax = int(n+0.5);
		}
    }

	printf("# Finish reading data...\n");

    /* calculate umb[][][] =  -beta*U */
	for(i = 0;i<nSimu;i++){
		for(j=0;j<nSimu;j++){	
	    	for(k=0;k<nFrame;k++){
	    		umb[i][j][k] = -beta*kappa[i]*0.5*( (nTw[j][k]-nStar[i]) * (nTw[j][k]-nStar[i])  ); 
	    	}
		}
    }

    /* Initialize the variables. */
	for (i = 0;i < nSimu-1; i++) {x[i] = 0.0;}

    /* Initialize the parameters for the L-BFGS optimization. */
    lbfgs_parameter_init(&param);
    param.linesearch = LBFGS_LINESEARCH_DEFAULT;

    /*
 *         Start the L-BFGS optimization; this will invoke the callback functions
 *                 evaluate() and progress() when necessary.
 *                      */
    
	printf("# L-BFGS begins...\n");
    ret = lbfgs(nSimu-1, x, &fx, evaluate, progress, umb, &param);

//	printf("# BFGS finish\n");
    /* Report the result. */
    printf("# L-BFGS optimization terminated with status code = %d\n", ret);
//    printf("# code = zero if the minimization process terminate without an error\n");
 /*
 * -1001 Rounding Error   =>  Normalize F and g in some way
 * 
 * */
//    printf("#  fx = %f\n",fx );
//    for(i = 0;i<nSimu-1;i++){
//		printf(" x[%d] = %.20f \n", i ,x[i]);
//    }

    /* Using x[i] i.e. delta ci to get ci*/

    //initial ci
	for(i = 0;i<nSimu;i++){ci[i] = 0.0;}    
    //ci = sum x from 0 to i-1
	for(i = 1;i<nSimu;i++){ci[i] = ci[i-1]+x[i-1];}
     
     /*using ci to get C*/
	for(i = 0;i<nSimu;i++){
		for(j = 0;j<nFrame;j++){
			for(k = 0;k<nSimu;k++){sumC[k] = log(nFrame)+umb[k][i][j]+ci[k];}		
			sumC_d[i+j*nSimu] = -logsumexp(nSimu,sumC);	
		}
    }
    C = -logsumexp(nSimu*nFrame,sumC_d);
	printf("# Normalized factor is C = %e\n",exp(C));

	sprintf(f,"WHAM/freeEnergy_window.dat");
    if(data ==NULL){printf("# error in creating data file %s\n# check if make the directory\n",f);}
	data=fopen(f,"w");
	for(i = 0;i<nSimu;i++){ fprintf(data,"c[%d] = %.4f ;\n",i,ci[i]);}
	fclose(data);

	//allocate memory for nTotal, pvN ,pvNt[]
    nTotal = (int *)malloc(sizeof(int)*(nMax+1));
    pvNt = (double *)malloc(sizeof(double)*(nTmax+1));
    pvN = (double *)malloc(sizeof(double)*(nMax+1));

	for(i = 0;i<(nMax+1);i++){
		pvN[i] = 0.0;
		nTotal[i] = 0;
	}
	for(i = 0;i<(nTmax+1);i++){ pvNt[i] = 0; }

    /*nTotal*/
    for(i=0;i<nSimu;i++){
        for(j = 0;j<nFrame;j++){
            nTotal[int( ni[i][j] + 0.5 )] += 1;
        }
    }
	sprintf(f,"WHAM/nTotal.dat");
	data = fopen(f,"w");
	if(data ==NULL){printf("# error in creating data file %s\n",f);}
	for(j=0;j<=nMax;j++){	fprintf(data,"%d\t%d\n",j,nTotal[j]); }
    fclose(data);
 
	/*pvN*/    
	for(i = 0;i<nSimu;i++){
		for(j = 0;j<nFrame;j++){
			for(k=0;k<nSimu;k++){sumC[k] = log(nFrame)+umb[k][i][j]+ci[k];}
			pvN[int(ni[i][j]+0.5)] += exp(C-logsumexp(nSimu,sumC));
		}
    }
    sprintf(f,"WHAM/pvN.dat");
    data = fopen(f,"w");
    if(data ==NULL){printf("# error in creating data file %s\n",f);}
	fprintf(data,"# n\tpvN\t-lnP\n");
    for(j=0;j<=nMax;j++){
        fprintf(data,"%d\t%e\t%e\n",j,pvN[j],-log(pvN[j]));
    }
    fclose(data);
    printf("# Pv(0) = %e; -logPv(0) = %e\n",pvN[0],-log(pvN[0]));
/*******************/

	/*pvNt*/    
	for(i = 0;i<nSimu;i++){
		for(j = 0;j<nFrame;j++){
			for(k=0;k<nSimu;k++){sumC[k] = log(nFrame)+umb[k][i][j]+ci[k];}
			pvNt[int(nTw[i][j]+0.5)] += exp(C-logsumexp(nSimu,sumC));
		}
	}  
    sprintf(f,"WHAM/pvNt.dat");
    data = fopen(f,"w");
	if(data ==NULL){printf("# error in creating data file %s\n",f);}
	fprintf(data,"# nt\tpvNt\t-lnP\n");
	for(j=0;j<nTmax;j++){
		fprintf(data,"%d\t%e\t%e\n",j,pvNt[j],-log(pvNt[j]));
	}
    fclose(data);

 //   printf("# finish pvN and pvNt");
  
    /*Calculate P_obs[window]*/
	double **pObs,*entropy;
	entropy = (double *)malloc(sizeof(double)*nSimu);
	pObs= (double **)malloc(sizeof(double *)*nSimu);
	for(i=0;i<nSimu;i++){
		pObs[i] = (double *)malloc(sizeof(double)*(nTmax+1));
	}
	for(i=0;i<nSimu;i++){
		for(j=0;j<=nTmax;j++){pObs[i][j]=0.0;}
	}

	for(i=0;i<nSimu;i++){
		for(j=0;j<nFrame;j++){
			pObs[i][int(nTw[i][j]+0.5)] += double(1.0/nFrame);
		}
    }

    double *sumpWHAM,**pWHAM;
	pWHAM = (double **)malloc(sizeof(double *)*nSimu);
	sumpWHAM = (double *)malloc(sizeof(double)*nSimu);
	for(i=0;i<nSimu;i++){
		sumpWHAM[i] = 0.0;
		pWHAM[i] = (double *)malloc(sizeof(double)*(nTmax+1));		
	}

	for(i=0;i<nSimu;i++){
		for(j=0;j<(nTmax+1);j++) pWHAM[i][j] = 0;
	}
    
    /*Calsulate S*/
    /*Normolize pWHAM by sum*/
	for(i=0;i<nSimu;i++){
		for(j=0;j<=nTmax;j++){
			if(pObs[i][j] > 0.0){
				pWHAM[i][j] = pvNt[j]*exp( - beta * kappa[i] * 0.5 * ( (j-nStar[i]) * (j-nStar[i]) ) + ci[i]); 
				sumpWHAM[i] += pWHAM[i][j];
			}
		}
    }
//    for(i=0;i<nSimu;i++){printf("# sumpWHAM[%d] = %f \n",i,sumpWHAM[i]);}

	for(i=0;i<nSimu;i++){
		entropy[i] = 0.0;
		for(j=0;j<=nTmax;j++){
			if(pObs[i][j] > 0.0){
				entropy[i] += pObs[i][j] * ( -log(pWHAM[i][j]) + log(sumpWHAM[i]) + log(pObs[i][j]) );
			}
		}
	}

	//display results
	sprintf(f,"WHAM/shanE.dat");
	data = fopen(f,"w");
	if(data ==NULL){printf("# error in creating data file %s\n",f);}
	for(j=0;j<nSimu;j++){
		fprintf(data,"%d\t%e\n",j,entropy[j]);
	}
    fclose(data);
    printf("# entropy = %e ; nTotal(0) = %d \n",entropy[0],nTotal[0]);

    /*display pObs (optional, turned off for now) */

/*
	for(i=0;i<nSimu;i++){
		sprintf(f,"WHAM/K%.3fN%d.dat",kappa[i],nStar[i]);	//need modify
		data = fopen(f,"w");
		if(data ==NULL){printf("# error in creating data file %s\n",f);}
		fprintf(data,"# N\tpObs\tPWHAM");
		for(j=0;j<=nTmax;j++){
			if(pObs[i][j]>0.0){
				fprintf(data,"%d\t%e\t%e\n",j,pObs[i][j],pWHAM[i][j]/sumpWHAM[i]);
			}
		}
		fclose(data);
    }
*/
	printf("# WHAM Finished...\n");

    return 0;
}

double logsumexp(int nx, void *y){
//calculate the log(sum exp(x)) of an double array x whose size is n
    int i;
    double max_x;
    double *x = (double *) y;
    double logsumexp,sumexp=0.0;
    //find the maximum of x
    max_x=x[0];
    for(i=1;i<nx;i++){
		if(x[i]>max_x){
	    	max_x = x[i];
		}
    }
    //sum_exp_(delta x)
    for(i=0;i<nx;i++){ 
		sumexp = sumexp + exp(x[i]-max_x);
    }
    //logsumexp = max + log_sum_exp(delta x)
    logsumexp = max_x + log(sumexp);

	return logsumexp;	//return logsumexp
}


//static double evaluate(
double evaluate(
    void *instance,
    const double *x,
    double *g,
    const int n,
    const double step
    )
{
	double *sumx;
	sumx = (double *)malloc(sizeof(double)*(n+1)); // g_i in Zhu, Hummer
	double fx = 0.0; // A^hat, the function to be minimized
	double fsum1 = 0.0; // negative of first term in eqn 22a
	double fsum2 = 0.0; // second term; -fsum1 + fsum2 = fx
	double *** umb = (double***) instance; // umb saves -\beta * U 
//	double **fln_max;
	int i,j,k,l,w = 0;
	double tem = log(nFrame);

	double *gsum2;
	gsum2 = (double *)malloc(sizeof(double)*(n+1)*nFrame); // 
	double *gsum2_d;
	gsum2_d = (double *)malloc(sizeof(double)*(n+1));
	double **gsum2_dlogsum;
	double *diffg;
	diffg = (double *)malloc(sizeof(double)*n);
	gsum2_dlogsum = (double **)malloc(sizeof(double *)*(n+1));
	for(i=0;i<=n;i++){gsum2_dlogsum[i] = (double *)malloc(sizeof(double)*nFrame);}
	
	sumx[0]=0.0;			//sumx[i] is c_i[i]
	for(i = 1;i<=n;i++){
		sumx[i] = sumx[i-1] +x[i-1];
	}	
	
	for(i = 1;i<=n;i++){
		fsum1 += sumx[i];
	}

	for(w=0;w<=n;w++){
	    for(l=0;l<nFrame;l++){
			for(k=0;k<=n;k++){
  		    	gsum2_d[k] = tem+umb[k][w][l]+sumx[k];
			}
			gsum2_dlogsum[w][l] = logsumexp((n+1),gsum2_d);
    	}
	}

	for(i = 0;i<=n;i++){
		for(j = 0;j<nFrame;j++){
			fsum2 += gsum2_dlogsum[i][j];
		} 
	}

	fx = -fsum1/(n+1)+fsum2/(nFrame*(n+1));

	//calculate g[n]
	//da/dg
	
	for(j = 1;j<=n;j++){
	    diffg[j] = 0;
	    for(w=0;w<=n;w++){
			for(l=0;l<nFrame;l++){
		    	gsum2[w*nFrame + l] = umb[j][w][l] - gsum2_dlogsum[w][l];
			}
	    }
	    diffg[j] = (exp(sumx[j]+logsumexp((n+1)*nFrame,gsum2))-1)/(n+1);
	}
	g[n-1] = diffg[n];
	for(i=n-2;i>=0;i--){
	    g[i] = g[i+1]+diffg[i+1];
	}
    return fx;
}

static int progress(
    void *instance,
    const double *x,
    const double *g,
    const double fx,
    const double xnorm,
    const double gnorm,
    const double step,
    int n,
    int k,
    int ls
    )
{
//    printf("#Iteration %d:\n", k);
//    printf("#  fx = %f\n",fx);
/*    for(int i =0;i<n;i++){
        printf(" x[%d] = %f\n",i, x[i]);
    }
    printf("#  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
    printf("\n");
*/    return 0;
}


