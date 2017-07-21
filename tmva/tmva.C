#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

using namespace std;

int main(int argc, char* argv[]){
//	char* variable_name = (argv[1]);	
//	char* type = (argv[2]);
	char cmd[50];
	sprintf(cmd,".x TMVAClassification_main.C(\"%s\",\"%s\",%i,%.2f,%i)",argv[1],argv[2],atoi(argv[3]),atof(argv[4]),atoi(argv[5]));
	cout<<cmd<<endl;
	gROOT->ProcessLine(cmd);
	return 0;
}
