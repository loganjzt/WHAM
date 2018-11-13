#include <cmath>
#include <cstdlib>
#include <iostream>
#include <cstring>
#include <fstream>
#include <sstream>
#include <ctime>
using namespace std;

int main () 
{
  int fnum;
  int nsims = 19;
  int nframes = 60001;
  int ncut = 10001;
  int ncount = nframes - ncut;
  double** N = (double**) malloc(nsims*sizeof(double*));
  for (int i=0; i<nsims; i++)
    N[i] = (double*) malloc(ncount*sizeof(double));
  double** Ntw = (double**) malloc(nsims*sizeof(double*));
  for (int i=0; i<nsims; i++)
    Ntw[i] = (double*) malloc(ncount*sizeof(double));
  double temp1, temp2;
  
  char format[] = "HC_K0.15_N%d.dat";
  char name[50];
  char str[100];
  FILE* infile;
  
  for (int i=0; i<nsims; i++)
  {
    fnum = -42 + i*14;
    sprintf (name, format, fnum);
    infile = fopen(name, "r");
    cout << name << endl;

    for (int j=0; j<4; j++)
      fgets (str, 100, infile);
    for (int k=0; k<ncut; k++)
      fscanf (infile, "%*f %*f %*f\n");
    for (int l=0; l<ncount; l++)
    {
      fscanf (infile, "%*f %lf %lf\n", &temp1, &temp2);
      N[i][l] = temp1;
      Ntw[i][l] = temp2;
    }
    fclose (infile);
  }

  // check inputs
  /*int n=6;
  for (int j=0; j<100; j++)
  {
    cout << N[n][j] << "\t" << Ntw[n][j] << endl;
  }*/

  int nbins = 300;
  double Nbin[nbins];
  double PN[nsims][nbins];
  double PNtot[nbins];
  int dN = 1;
  int bin;
  for (int i=0; i<nsims; i++)
  {
    for (int k=0; k<nbins; k++)
    {
      Nbin[k] = 0.0;
      PN[i][k] = 0.0;
    }
  }

  for (int i=0; i<nsims; i++)
  {
    for (int j=0; j<ncount; j++)
    {
      bin = int(N[i][j]/dN);
      PN[i][bin] += 1.0;
      PNtot[bin] += 1.0;
    }
  }
    
  for (int k=0; k<nbins; k++)
  {
    Nbin[k] = k*dN;
    cout << Nbin[k] << "\t" << PNtot[k] << "\t";
    for (int i=0; i<nsims; i++)
    {
      cout << PN[i][k] << "\t";
    }
    cout << "\n";
  }
    
  return 0;
}
 
