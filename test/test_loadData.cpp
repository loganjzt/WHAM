#include <fstream>	// std::ifstream
#include <cstdio>	// IO operations , ig. printf
#include <cmath>	// log functions
#include "/home/zhitongj/Git/wham/misc_fn.cpp"


	using namespace std;


// test loadData fn

	double *logebfk,*difflogebfk;
	double *mu,*kappa,*x1Star;

	int **x1;
	double **x1Tilda;

	int *nSkip, *nData;	
	double *lognData;
	int nSimu,nDataMax;
	double temperature;
	double kBT;
	char **fnarray;

#include "/home/zhitongj/Git/wham/loadData.cpp"

int main(int argc,char* argv[]){

	loadData(argv[1]);
	return 0;
}



