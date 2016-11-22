#include "TFile.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"
#include "CreateTree_tmva_all.h"
#include "CreateTree_tmva_all.C"
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <string>
#include <stdlib.h>

using namespace std;

int main(int argc, char* argv[]){
		TString input_filename = std::string(argv[1]);
		TString file_tag = std::string(argv[2]);
		TString region  = std::string(argv[3]);
		
		TFile *f = TFile::Open(input_filename);
		if (f!=0){	
			CreateTree_tmva_all	*c = new CreateTree_tmva_all(0,input_filename);
			c->Loop(file_tag, region);
		}
	return 0;
}
