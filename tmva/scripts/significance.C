#include <TH1D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TLorentzVector.h"
#include <iomanip>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "TFile.h"
#include "TROOT.h"

//int main(){
void significance(){
	gROOT->ProcessLine(".x /afs/cern.ch/work/n/nchernya/setTDRStyle.C");

	Double_t hist_integrals;
	TString file_names = "/shome/nchernya/Hbb/tmva/main_mva/output/NmN/v14/SignN_1/TMVA_main_v14_Data_Nm1_all_double_csv_jet5.root";

		TFile *file = new TFile(file_names);
		file->cd("Method_BDT/BDTG");
	//	file->ls();
	//	TH1D *hist_S = (TH1D*) MVA_BDTG_S->Clone();
		hist_S->Scale(1./hist_S->Integral()*(132.8-52.3)); 
		hist_S->Scale(1./hist_S->Integral()*(75.3)); 
		TH1D *hist_B = (TH1D*) MVA_BDTG_B->Clone(); 
		hist_B->Scale(1./hist_B->Integral()*(2040650.-119955.)); 
	//	hist_B->Scale(1./hist_B->Integral()*(582247.)); 
		
		double END = hist_B->GetBinCenter(hist_B->FindLastBinAbove(0.)); //right end of BDT distibution
		double bin=0.;
		double s1, b1;
		int i=0;
		double max=0;	
		do	{
			s1=hist_S->GetBinContent(i+1);
			b1=hist_B->GetBinContent(i+1);
			bin=(double) hist_S->GetBinCenter(i+1+1);
			if (b1!=0) max += pow(s1,2)/b1;
			i++;
		} while (bin < END);
	
	ofstream out;
	out.open("double_csv_jet5.txt");
	out<<max<<endl;
	cout<<max<<endl;

}
