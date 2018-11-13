//
// test call_bfgs.cpp
//
#include "/home/pamish/packages/utils/lbfgs_install/include/lbfgs.h"
#include <stdio.h>

	using namespace std;

// test call_bfgs fn

	double *logebfk,*difflogebfk;
	double ***logebWk;
	int nSimu;
	int *nData;
	double *lognData;
	int nDataMax;	

#include "/home/zhitongj/Git/wham/call_bfgs.cpp"

int main(int argc,char* argv[]){

	int call_bfgs();
	call_bfgs();

	return 0;
}


