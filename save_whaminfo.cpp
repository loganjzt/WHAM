//
//	save_whaminfo.cpp
//
#include <fstream>
#include <stdio.h>
#include <iostream>
#include <string>
#include <cstring>
#include <stdlib.h>
#include <cstdio> // fopen

	using namespace std;

int save_whaminfo(char* infoFileName){
  
	FILE *infoFileO;	

	printf("# save wham info to %s\n",infoFileName);   

	infoFileO = fopen(infoFileName,"w");
	if(infoFileO == NULL){
		printf("# Error : unable to create save whaminfo file\n ");
		return 1;
	}
   
	fprintf (infoFileO, "#nsim        temperature \n");
	fprintf (infoFileO, "%d\t%f\n",nSimu,temperature);
	fprintf (infoFileO, "# srno   nskip   ndata  mu   Kappa  Nstar  logebfk filename \n");

	for(int i = 0; i < nSimu; i++) {
		fprintf (infoFileO, "%d\t%d\t%d\t%f\t%f\t%f\t%18f\t%s\n",i+1,nSkip[i],nData[i],kBT*mu[i],kBT*kappa[i],x1Star[i],logebfk[i],fnarray[i]);
	}
  
	fclose (infoFileO);

	return 0;
}


