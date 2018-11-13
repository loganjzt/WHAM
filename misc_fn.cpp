//
//  misc_fn.c
//  
//
//  Created by Amish J Patel on 4/23/13.
//
//  Modified by Erte Xi on 5/19/2015 
//
//  Modified by Zhitong Jiang on 11/09/2018

#include <stdio.h>
#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <cmath>
#include <math.h>
#include <stdlib.h>

#ifndef MISC_FN
#define MISC_FN

// Potentials
// Move this part to data2ebWk fn
//double linear_umbrella(double n, double m, double k, double nc)
//{
//   return(  m * n  +  0.5 * k * (n - nc) * (n - nc)  );
//}


// LOGSUMEXP
double logsumexp(int nx, double *x)
{
   double sumexp, maxn;
   double max_double_1darray(int,double*);
   
   maxn = max_double_1darray(nx,x);
   sumexp = 0.0;
   for(int i=0;i<nx;i++){
      if ( ( maxn - x[i] ) < 600.0 ) sumexp = sumexp + exp( x[i] - maxn );
   }
   return( log(sumexp) + maxn );
}

// MAXMIN
double max_double_1darray(int nx, double *x)
{
   double maxx=x[0];
   for(int i=1;i<nx;i++) {
      if(x[i]>maxx) maxx=x[i];
   }
   return (maxx);
}

double min_double_2darray(int nx, int ny, double **x)
{
   double minx=x[0][0];
   for(int i=0;i<nx;i++) {
      for(int j=0;j<ny;j++) {
         if(x[i][j]<minx) {
            minx=x[i][j];
         }
      }
   }
   return(minx);
}

double max_double_2darray(int nx, int ny, double **x)
{
   double maxx=x[0][0];
   for(int i=0;i<nx;i++) {
      for(int j=0;j<ny;j++) {
         if(x[i][j]>maxx) {
            maxx=x[i][j];
         }
      }
   }
   return(maxx);
}

int max_int_1darray(int nx, int *x)
{
   int maxx=x[0];
   for(int i=1;i<nx;i++) {
      if(x[i]>maxx) maxx=x[i];
   }
   return (maxx);
}

int min_int_2darray(int nx, int ny, int **x)
{
   int minx=x[0][0];
   for(int i=0;i<nx;i++) {
      for(int j=0;j<ny;j++) {
         if(x[i][j]<minx) {
            minx=x[i][j];
         }
      }
   }
   return(minx);
}

int max_int_2darray(int nx, int ny, int **x)
{
   int maxx=x[0][0];
   for(int i=0;i<nx;i++) {
      for(int j=0;j<ny;j++) {
         if(x[i][j]>maxx) {
            maxx=x[i][j];
         }
      }
   }
   return(maxx);
}


// CREATING ARRAYS FOR GLOBAL USE

int** Create2dIntArray(int x, int y)
{
	int** temp;
	temp = new int* [x];
	for (int i=0; i < x ; i++)
        {
	  temp[i] = new int[y];
          for (int j=0; j<y; j++)
          {
            temp[i][j] = 0;
          }
	}
        return temp;
}

double** Create2dDoubleArray(int x, int y)
{
        double** temp;
        temp = new double*[x];
        for(int i=0;i<x;i++)
        {
          temp[i] = new double[y];
          for (int j=0; j<y; j++)
          {
            temp[i][j] = 0.0;
          }
        }
        return temp;
}

char** Create2dCharArray(int r, int c)
{
	char** temp;
	temp = new char*[r];
	for(int x=0;x<r;x++)
		temp[x] = new char[c];
	return temp;
}

double*** Create3dDoubleArray(int x, int y, int z)
{
	double*** temp;
	temp = new double**[x];
	for (int i=0; i<x; i++) 
        {
          temp[i] = new double*[y];
          for (int j=0; j<y; j++) 
          {
            temp[i][j] = new double[z];
            for (int k=0; k<z; k++)
            {
              temp[i][j][k] = 0.0;
            }
          }
        }
	return temp;
}

// Free the memory of the created arrays
void Free2dDoubleArray(double** the_array, int x, int y)
{
  for ( int i=0; i<x; i++ )
  {
    delete [] the_array[i];
  }
  delete [] the_array;
}

#endif

