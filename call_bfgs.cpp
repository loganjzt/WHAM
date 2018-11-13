//
// call_bfgs.cpp
//
//	Modified from AJP by Zhitong Nov 12 2018
//
#include "/home/pamish/packages/utils/lbfgs_install/include/lbfgs.h"
#include <stdio.h>
#include "/home/zhitongj/Git/wham/misc_fn.cpp"

	using namespace std;

int call_bfgs(){
	
	double evaluate(void *instance,const double *x,double *g,const int n,const double step);
	int progress(void *instance,const double *x,const double *g,const double fx,const double xnorm,const double gnorm,const double step,int n,int k,int ls);

	const int nDim = nSimu - 1;
	int ret = 0;
	double fx;
	double *x = lbfgs_malloc(nDim);
	lbfgs_parameter_t param;

	if( x == NULL){
		printf("# Error: Failed to allocate memory for lbfgs variables.\n");
		return 1;
	}

	// Initialize lbfgs variables
	for(int i = 0; i < nDim ; i++) x[i] = difflogebfk[i];
	lbfgs_parameter_init(&param);
	param.linesearch = LBFGS_LINESEARCH_DEFAULT;

	printf("# L_BFGS begin...\n");
	// Start L-BFGS optimization;
    ret = lbfgs(nDim, x, &fx, evaluate, progress, logebWk, &param);

	// Report the result
	printf("# L-BFGS optimization terminated with status code = %d\n", ret);
	// 0 is success, 2 LBFGS_ALREADY_MINIMIZED
	// errors: -1022 LBFGSERR_OUTOFMEMORY, -1001 LBFGSERR_ROUNDING_ERROR,
	// -1003 LBFGSERR_OUTOFINTERVAL,-998 LBFGSERR_MAXIMUMLINESEARCH,
	// -997 LBFGSERR_MAXIMUMITERATION, -996 LBFGSERR_WIDTHTOOSMALL, -994 LBFGSERR_INCREASEGRADIENT,
	
   	for (int i = 0; i < nDim; i++) printf("# x[%d] = %g \n",i+1,x[i]);
	printf("# fx=%g\n",fx);

	for(int i = 0 ; i < nDim ; i++) difflogebfk[i] = x[i];

	for(int i = 0; i < nSimu; i++) {
		logebfk[i] = 0.0;
		for(int j = 0; j < i ; j ++) logebfk[i] += difflogebfk[j];
	}

	lbfgs_free(x);

	return 0;
}


double evaluate(
    void *instance,
    const double *x,
    double *g,
    const int n,
    const double step
    )
{
	double *** umb = (double***) instance; // umb saves -\beta * U 
	double *sumx = new double [n+1];	// g_i in Zhu, Hummer, eq 21b

	double fx = 0.0; // A^hat, the function to be minimized
	double fsum1 = 0.0; // negative of first term in eqn 22a
	double fsum2 = 0.0; // second term; -fsum1 + fsum2 = fx

	int i,j,k,l,w = 0;

	double *gsum2 = new double [nDataMax*(n+1)];
	//gsum2 = (double *)malloc(sizeof(double)*(n+1)*nFrame);
	double *gsum2_d = new double [n+1];
	//gsum2_d = (double *)malloc(sizeof(double)*(n+1));
	double *diffg = new double [n];

	double **gsum2_dlogsum = new double* [n+1];
	gsum2_dlogsum = (double **)malloc(sizeof(double *)*(n+1));
	for(i = 0 ; i < n + 1 ; i++){
		//gsum2_dlogsum[i] = (double *)malloc(sizeof(double)*nFrame);
		gsum2_dlogsum[i] = new double [nDataMax];
	}
	
	sumx[0]=0.0;			//sumx[i] is c_i[i]
	for(i = 1; i <= n; i++){
		sumx[i] = sumx[i-1] +x[i-1];
	}	
	
	for(i = 1; i <= n; i++){
		fsum1 += sumx[i] * nData[i] ;
	}

	for(w = 0; w <= n; w++){
	    for(l=0 ; l < nData[w] ;l++){
			for(k = 0 ; k <= n; k++){
  		    	gsum2_d[k] = lognData[k]+umb[k][w][l]+sumx[k];
			}
			gsum2_dlogsum[w][l] = logsumexp((n+1),gsum2_d);
    	}
	}

	for(i = 0; i<= n; i++){
		for(j = 0; j < nData[i] ; j++){
			fsum2 += gsum2_dlogsum[i][j];
		} 
	}

	fx = -fsum1 / ( nDataMax*(n+1) ) + fsum2/(nDataMax*(n+1));

	//calculate g[n]
	// dA/dg
	int ii = 0;	
	for(j = 1; j <= n; j++){
	    diffg[j] = 0;
		ii = 0;
	    for(w = 0; w <= n; w++){
			for(l = 0; l < nData[j]; l++){
		    	gsum2[ii] = umb[j][w][l] - gsum2_dlogsum[w][l];
				ii++;
			}
	    }
	    diffg[j] = ( exp( sumx[j] + logsumexp(ii,gsum2) ) -1.0 ) / (n+1);
	}

	g[n-1] = diffg[n];
	for(i=n-2;i>=0;i--){
	    g[i] = g[i+1]+diffg[i+1];
	}
	for(i = 0; i < n; i++) g[i] = g[i] ;

	// release memory
	delete [] sumx;
	delete [] diffg;
	delete [] gsum2;
	delete [] gsum2_d;
	for(int i =0; i < n+1 ; i++){
		delete [] gsum2_dlogsum[i];
	}

    return fx;
}

int progress(
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
	printf("#Iteration %d:\n", k);
	printf("#  fx = %f\n",fx);
/*
	for(int i =0;i<n;i++){
        printf(" x[%d] = %f\n",i, x[i]);
    }
	printf("#  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
	printf("\n");
*/    
	return 0;
}

