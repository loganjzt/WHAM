#define LBFGS_FLOAT 64
#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <cmath>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "../uwham/wham.h"
using namespace std;

// ALWAYS CHECK UNITS OF KBT IN BOTH LOGWHAM AND GETPVN OR POT_OF_MF AS WELL AS THE INPUT FORMAT FOR AND ORDER OF NREG AND NTILDA !!!!

int nsim, ndatamax, ndatatot, ijoint;
double temperature, kBT;
int *nskip, *ndata; // [nsim]
double *mu, *Nstar, *Kappa, *logebfk, *difflogebfk, *logndata; // [nsim]
int **Nreg; // dims [nsim][ndatamax]
double **Ntilda; // dims [nsim][ndatamax]
double ***logebWk; //dims [nsim][nsim][ndatamax] :  ebWk(k,i,l) will store exp(-beta*Wk(Ntilda(i,l))
char** fnarray;

int Nregmin, Nregmax, nNreg;
int *Nregval, *nNregpts; //dims [nNreg];
double *PvNreg, *logPvNreg; // dims [nNreg];
double dNtilda;
int Ntildamin, Ntildamax, nNtilda;
int *nNtildapts; //dims [nNtilda];
double *Ntildaval, *PvNtilda, *logPvNtilda; // dims [nNtilda];

int** nNregNtildapts; // dims [nNreg][nNtilda]
double **PvNregNtilda, **logPvNregNtilda; // dims [nNreg][nNtilda]


int main ()
{
 
  void load_ebWk( string,int,double(*)(double,double,double,double),double,int );
//Usage:load_ebWk(file_whaminfo,icol,potentialfn,dNtilda,ijoint)
// icol is the column number with Nreg data, icol+1 has the Ntilda data
   int call_bfgs(int);
   void save_whaminfo(string);
   void get_PvN(string,double,int);
//Usage:get_PvN(file_whaminfo,Ntildafine,isave)
   void check_PvN(string,double(*)(double,double,double,double),int);
//Usage: void check_PvN(file_whaminfo,potentialfn,isave)
  double linear_umbrella(double,double,double,double);
   
   int ijoint=0;
   load_ebWk("Ndata",2,linear_umbrella,1.0,ijoint);
   
   call_bfgs(0);
   
   int isave=1;
   if(isave==1) save_whaminfo("Ndata");
   
   get_PvN("Ndata",1.0,isave);

   check_PvN("Ndata",linear_umbrella,isave);
   
  return 0;
}

