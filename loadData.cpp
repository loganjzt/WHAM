#include <fstream>	// std::ifstream
#include <cstdio>	// IO operations , ig. printf
#include <cmath>	// log functions
#include "/home/zhitongj/Git/wham/misc_fn.cpp"

	using namespace std;


// test loadData fn

int loadData(char* infoFileName)
{	// Modified from AJP's WHAM code

/*
	double max_double_1darray(int, double*);
	int** Create2dIntArray(int,int);
	double** Create2dDoubleArray(int,int);	
	char** Create2dCharArray(int,int);
	double*** Create3dDoubleArray(int,int,int);
	int min_int_2darray(int, int, int**);
	int max_int_2darray(int, int, int**);
	double min_double_2darray(int, int, double**);
	double max_double_2darray(int, int, double**);
  */ 
	int windowN;	//	# of window
	double time;
	double x1tmp;

	// read WHAMINFO file
	ifstream fin;
	fin.open(infoFileName);
	fin.ignore(1000,'\n'); 
	fin >> nSimu;
	fin.ignore(1000,'\t'); 	
	fin >> temperature;
	fin.ignore(1000,'\n'); 
	fin.ignore(1000,'\n'); 
	printf("# %d windows, temperature = %.3f\n",nSimu,temperature);

	// nSimu dependent memory allocation
	nSkip = new int[nSimu]; // number of data points to be skipped
	nData = new int[nSimu]; // number of data points to be analyzed
	lognData = new double[nSimu];	
	
	mu = new double[nSimu];
	kappa = new double[nSimu];
	x1Star = new double[nSimu];

	logebfk = new double[nSimu];
	difflogebfk = new double[nSimu-1];
 
	// char (*fnarray)[500] = new char[nSimu][500];
	fnarray = Create2dCharArray(nSimu,500);
	// resume reading WHAMINFO file
	for(int i=0;i<nSimu;i++) {
		fin >> windowN;
		fin.ignore(1000,'\t');
		fin >> nSkip[i];
		fin.ignore(1000,'\t');
		fin >> nData[i];
		lognData[i] = log( nData[i] );
		fin.ignore(1000,'\t');
		fin >> mu[i];
		fin.ignore(1000,'\t'); 
		fin >> kappa[i];
		fin.ignore(1000,'\t'); 
		fin >> x1Star[i];
		fin.ignore(1000,'\t'); 
		fin >> logebfk[i];
		fin.ignore(1000,'\t'); 
		fin >> fnarray[i];
		fin.ignore(1000,'\n'); 
		//cout << windowN << "\t" << nSkip[i] << "\t" << nData[i] << "\t" << mu[i] << "\t" << kappa[i] << "\t" << x1Star[i] << "\t" << logebfk[i] << "\t" << fnarray[i] << endl;
		printf("%d\t%d\t%d\t%f\t%f\t%f\t%f\t%s\n",windowN,nSkip[i],nData[i],mu[i],kappa[i],x1Star[i],logebfk[i],fnarray[i]);
	}
	fin.close();

	// check input
	if(windowN!=nSimu){
		printf("# Info file read incomplete!\n# Aborting...\n");
		abort();
	}	
	// finish reading WHAMINFO file

	// change all unit to kT 
	kBT = 8.314 * temperature / 1000.0;
	for(int i = 0 ; i < nSimu ; i++ ){
		kappa[i] /= kBT;
		mu[i] /= kBT;
	}

	// get difflogebfk from logebfk
	//  use the first window as ref state
	for(int i = nSimu ; i > 0 ; i-- ) logebfk[i-1] = logebfk[i-1] - logebfk[0];
	for(int i = 0 ; i < nSimu - 1 ; i++ ) difflogebfk[i] = logebfk[i+1] - logebfk[i];
	

	// allocate memory for nSimu and nData dependent 
	nDataMax = max_int_1darray(nSimu,nData);

	x1 = Create2dIntArray( nSimu , nDataMax );
	x1Tilda = Create2dDoubleArray( nSimu , nDataMax );

	// read x1 and x1Tilda from each windows
	for(int i = 0; i < nSimu; i ++ ){
		ifstream fin;
		fin.open(fnarray[i]);	// open a data file
		for(int j = 0; j < nSkip[i]; j ++) fin.ignore(10000,'\n');	 // skip some data from eq or block average
		for(int j = 0; j < nData[i]; j ++){
			fin >> time ; 
			fin.ignore(100,'\t');
			fin >> x1tmp;
			x1[i][j] = int(x1tmp+0.5);
			fin.ignore(100,'\t');
			fin >> x1Tilda[i][j];
			fin.ignore(1000,'\n');
			if(j == 0 ){
				printf("# %d\t%g\t%g\t%d\n",i + 1  , x1Star[i], x1Tilda[i][j] , x1[i][j] );
			}
		}
		fin.close();
	}

	// Finish reading x1 and x1Tilda data

	return 0;
}


