#include <stdio.h>
#include <iomanip>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2.h"
#include "TTree.h"
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TMarker.h>
#include <THStack.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TCut.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TF1.h>
#include "TROOT.h"
#include "TPaveText.h"
#include "TPad.h"
#include "TGaxis.h"
#include <boost/format.hpp>
#include <array>
#include <algorithm>
#include <iterator>

#define SWAP(A, B) { Float_t t = A; A = B; B = t; }
#define SWAP2(A, B) { B.swap(A); }

void bubblesort(Float_t *a, std::string *str, int n)
{
  int i, j;
  for (i = n - 1; i >= 0; i--)
  {
    for (j = 0; j < i; j++)
    {
      if (a[j] < a[j + 1]){
        SWAP( a[j], a[j + 1] );
	//	cout<<str[j]<<"   "<<str[j+1]<<endl;
		SWAP2( str[j], str[j+1]);
	//	cout<<str[j]<<"   "<<str[j+1]<<endl;
		}
    }
  }
}


Float_t Discr(TH1F *h1, TH1F *h2){  
  Int_t h1s = h1->GetNbinsX();
  Int_t h2s = h2->GetNbinsX();
  if (h1s != h2s) {
  	printf("h1 size %d   h2 size %d \n",h1s,h2s);
  	return -1;
  }
  if (h1->Integral()!=0 ) h1->Scale(1./h1->Integral());
  if (h2->Integral()!=0 ) h2->Scale(1./h2->Integral());
  Float_t adiff = 0;
  for (Int_t i=0;i<h1s;i++){
	adiff +=  TMath::Abs(h1->GetBinContent(i) - h2->GetBinContent(i));
	//cout << "h1 bin " << i << " " << h1->GetBinContent(i) << endl;
	}
  return adiff/2;
}

using namespace std;

int main(int argc, char* argv[]){
//int analyzer_stackRatio(){
gROOT->ProcessLine(".x /afs/cern.ch/work/n/nchernya/setTDRStyle.C");
int set_type=atoi(argv[1]); // 0 - analysis, 1 - control region , top
 
const int nfiles  = 5; //9;  //11;
//TString leg_names[nfiles] = {"Data","VBF Z #rightarrow ll (#times 10)","WW + jets","ZZ + jets","WZ + jets", "t#bar{t}", "Z + jets"};
//TString leg_names[nfiles] = {"Data","tZq #rightarrow ll","TTZ #rightarrow ll #nu#nu","W(l#nu) + jets","WW + jets","ZZ + jets","WZ + jets","Single Top", "t#bar{t}", "Z + jets","VBF Z #rightarrow ll"};


//TString leg_names[nfiles] = {"Data","W(l#nu) + jets","WW + jets","ZZ + jets","WZ + jets","Single Top", "t#bar{t}", "Z + jets (MDG)","VBF Z #rightarrow ll"};
//TString leg_names[nfiles] = {"Data","W(l#nu) + jets","WW + jets","ZZ + jets","WZ + jets","Single Top", "t#bar{t}", "Z + jets (AMC)","VBF Z #rightarrow ll"};
TString leg_names[nfiles] =  {"Data","VV","Top quark", "Z + jets","EW Zjj"};
TString leg_names_pas[nfiles] = {"Data","VV","Top quark", "Z + jets","EW Zjj"};
//TString leg_names[nfiles] = {"Data","WW + jets","ZZ + jets","WZ + jets","Single Top", "t#bar{t}", "Z + jets (AMC)","VBF Z #rightarrow ll"};
TString set_names[2] = {"Dimuon","Dielectron"}; 
//if (set_type==0)leg_names[0] = "Data (DoubleB)";
//if (set_type==1)leg_names[0] = "Data (SingleB)";

/*
TString file_names[nfiles] = {"SingleMuon","EWK_LLJJ","WW","ZZ","WZ","TT","DYJetstoLL_HT"};
TString file_names_QCDup[nfiles] = {"SingleMuon","EWK_LLJJ","WW","ZZ","WZ","TT","DYJetstoLL_HT"};
TString file_names_QCDdown[nfiles] = {"SingleMuon","EWK_LLJJ","WW","ZZ","WZ","TT","DYJetstoLL_HT"};
TString file_names_JESup[nfiles] = {"SingleMuon","EWK_LLJJ","WW","ZZ","WZ","TT","DYJetstoLL_HT"};
TString file_names_JESdown[nfiles] = {"SingleMuon","EWK_LLJJ","WW","ZZ","WZ","TT","DYJetstoLL_HT"};

TString file_names[nfiles] = {"SingleMuon","tZq_ll","TTZToLLNuNu","WJetsToLNu","WW","ZZ","WZ","ST","TT","DYJetstoLL_HT","EWK_LLJJ"};
TString file_names_QCDup[nfiles] = {"SingleMuon","tZq_ll","TTZToLLNuNu","WJetsToLNu","WW","ZZ","WZ","ST","TT","DYJetstoLL_HT","EWK_LLJJ"};
TString file_names_QCDdown[nfiles] = {"SingleMuon","tZq_ll","TTZToLLNuNu","WJetsToLNu","WW","ZZ","WZ","ST","TT","DYJetstoLL_HT","EWK_LLJJ"};
TString file_names_JESup[nfiles] = {"SingleMuon","tZq_ll","TTZToLLNuNu","WJetsToLNu","WW","ZZ","WZ","ST","TT","DYJetstoLL_HT","EWK_LLJJ"};
TString file_names_JESdown[nfiles] = {"SingleMuon","tZq_ll","TTZToLLNuNu","WJetsToLNu","WW","ZZ","WZ","ST","TT","DYJetstoLL_HT","EWK_LLJJ"};

*/
/*
TString file_names[nfiles] = {"SingleMuon","WJetsToLNu","WW","ZZ","WZ","ST","TT","DYJetstoLL_amc","EWK_LLJJ"};
TString file_names_QCDup[nfiles] = {"SingleMuon","WJetsToLNu","WW","ZZ","WZ","ST","TT","DYJetstoLL_amc","EWK_LLJJ"};
TString file_names_QCDdown[nfiles] = {"SingleMuon","WJetsToLNu","WW","ZZ","WZ","ST","TT","DYJetstoLL_amc","EWK_LLJJ"};
TString file_names_JESup[nfiles] = {"SingleMuon","WJetsToLNu","WW","ZZ","WZ","ST","TT","DYJetstoLL_amc","EWK_LLJJ"};
TString file_names_JESdown[nfiles] = {"SingleMuon","WJetsToLNu","WW","ZZ","WZ","ST","TT","DYJetstoLL_amc","EWK_LLJJ"};*/

TString file_hist_names[nfiles] = {"data","VV","top","DY","EWK_LL"};
TString file_names[nfiles]={};

int bg_begin;
int qcd_begin=18;
 //bg_begin=2;
 bg_begin=1;
/*
int FILLCOLOR[nfiles] = {1,kRed-4,kSpring+8, kSpring+5, kGreen-3, kBlue-4,kOrange-2};
int LINECOLOR[nfiles] = {1,kRed-4,kSpring+8, kSpring+5, kGreen-3,kBlue-4,kOrange-2};
int LINESTYLE[nfiles] = {1,1,1,1,1,1,1};
int LINEWIDTH[nfiles] = {1,3,1,1,1,1,1};
int FILLSTYLE[nfiles] = {1001,1001,1001,1001,1001,1001,1001};

int FILLCOLOR[nfiles] = {1,kTeal+7,kTeal-1,kSpring+7,kSpring+8, kSpring+5, kGreen-3, kAzure+1,kBlue-4,kOrange-2,kRed-4};
int LINECOLOR[nfiles] = {1,kTeal+7,kTeal-1,kSpring+7,kSpring+8, kSpring+5, kGreen-3,kAzure+1,kBlue-4,kOrange-2,kRed-4};
int LINESTYLE[nfiles] = {1,1,1,1,1,1,1,1,1,1,1};
int LINEWIDTH[nfiles] = {1,1,1,1,1,1,1,1,1,1,3};
int FILLSTYLE[nfiles] = {1001,1001,1001,1001,1001,1001,1001,1001,1001,1001,1001};*/
//int FILLCOLOR[nfiles] = {1,kSpring+7,kSpring+8, kSpring+5, kGreen-3, kAzure+1,kBlue-4,kOrange-2,kRed-4};
//int LINECOLOR[nfiles] = {1,kSpring+7,kSpring+8, kSpring+5, kGreen-3,kAzure+1,kBlue-4,kOrange-2,kRed-4};
int FILLCOLOR[nfiles] = {1,kGreen-3, kAzure+1,kOrange-2,kRed-4};
int LINECOLOR[nfiles] = {1,kGreen-3,kAzure+1,kOrange-2,kRed-4};
int LINESTYLE[nfiles] = {1,1,1,1,1};
int LINEWIDTH[nfiles] = {1,1,1,1,3};
int FILLSTYLE[nfiles] = {1001,1001,1001,1001,1001};
/*
int FILLCOLOR[nfiles] = {1,kSpring+8, kSpring+5, kGreen-3, kAzure+1,kBlue-4,kOrange-2,kRed-4};
int LINECOLOR[nfiles] = {1,kSpring+8, kSpring+5, kGreen-3,kAzure+1,kBlue-4,kOrange-2,kRed-4};
int LINESTYLE[nfiles] = {1,1,1,1,1,1,1,1};
int LINEWIDTH[nfiles] = {1,1,1,1,1,1,1,3};
int FILLSTYLE[nfiles] = {1001,1001,1001,1001,1001,1001,1001,1001};
*/
	
//int order[nfiles] = {0,1,2,3,4,5,6,7,8,9,10};// number in file_names array as they should appear in the legend
int order[nfiles] = {0,1,2,3,4};// number in file_names array as they should appear in the legend
int order_legend[nfiles]; 
for (int i=0;i<nfiles;i++){
	order_legend[order[i]]=i;
}
	
TString data_name[2] = {"SingleMuon","SingleElectron"};

TString set[3]={"_mu","_el"}; 

for (int i=0;i<nfiles;i++){
	if (i==0) {
		if (set_type==0)	file_names[i] = ("/afs/cern.ch/work/n/nchernya/VBFZll/forZpt_plot_root/new_files/histos_EWK_LL_ptZ_mu_SMfromATGC_dataOnly_PoissonErrors.root");
		if (set_type==1)	file_names[i] = ("/afs/cern.ch/work/n/nchernya/VBFZll/forZpt_plot_root/new_files/histos_EWK_LL_ptZ_el_SMfromATGC_dataOnly_PoissonErrors.root");
	} else{
  if (set_type==0) file_names[i] = ("/afs/cern.ch/work/n/nchernya/VBFZll/forZpt_plot_root//histos_EWK_LL_ptZ_mu_SMfromATGC_added.root");
  if (set_type==1) file_names[i] = ("/afs/cern.ch/work/n/nchernya/VBFZll/forZpt_plot_root/histos_EWK_LL_ptZ_el_15bins_added.root");
//  if (set_type==0) file_names[i] = ("/afs/cern.ch/work/n/nchernya/VBFZll/forZpt_plot_root/new_files/histos_EWK_LL_ptZ_mu_SMfromATGC_BDTcut_added.root");
  //if (set_type==1) file_names[i] = ("/afs/cern.ch/work/n/nchernya/VBFZll/forZpt_plot_root/new_files/histos_EWK_LL_ptZ_el_SMfromATGC_BDTcut_added.root");
  }
 //  file_names[i].Prepend("root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat//store/user/nchernya/VBFZll/plotterOutput/v25/");
//	else  file_names[i].Prepend("/afs/cern.ch/work/n/nchernya/VBFZll/plotter/output_root/");
//	file_names[i].Append(set[set_type]);
//	file_names_QCDup[i] = file_names[i];
//	file_names_QCDdown[i] = file_names[i];
//	file_names_JESup[i] = file_names[i];
//	file_names_JESdown[i] = file_names[i];
//	}
//	file_names[i].Append(".root");
//	file_names_QCDup[i].Append(".root");
//	file_names_QCDdown[i].Append(".root");
//	file_names_JESup[i].Append(".root");
//	file_names_JESdown[i].Append(".root");
}
TString trigger[2] = {"", ""};
//TString dir_name= "plots_mdg_ht_bdt_alldata_v25_reaod";
//TString dir_name= "plots_amc_bdt_alldata_v25_reaod";
//TString dir_name= "plots_amc_bdt_axis2_v25_reaod";
//TString dir_name= "plots_amc_bdt_alldata4qglnorm_v25_reaod_pas_herwig";
//TString dir_name= "paper/plots_herwig";
TString dir_name= "plots_29112017";
//TString dir_name= "pas_dir/plots_amc_bdt_alldata4qglnorm_v25_reaod_pas";
///////TString dir_name= "test_qgl_look";
//TString dir_name= "plots_mdg_ht_bdt_axis2_v25_reaod";
//TString dir_name= "plots_amc_bdt_v25";
//TString dir_name= "plots_amc_herwig_bdt_v25";
dir_name.Append(set[set_type]+"/");
//TString dir_name = "plots_amc/";
Float_t lumi = 35900;

	


//TLegend *leg = new TLegend(0.77,0.65,0.92,0.9); //without writing about SF   for AN
//////TLegend *leg = new TLegend(0.73,0.60,0.92,0.85); //without writing about SF   for PAS
TLegend *leg = new TLegend(0.7,0.50,0.92,0.85); //without writing about SF   for PAS
//TLegend *leg = new TLegend(0.73,0.50,0.92,0.85); //without writing about SF   for final paper but cw, noot cw/l2
//TLegend *leg = new TLegend(0.7,0.60,0.92,0.85); //without writing about SF   for PAS with linear qgl
leg->SetFillColor(0);
leg->SetBorderSize(0);
leg->SetTextFont(42);
leg->SetTextSize(0.025);
TLegend *leg2 = new TLegend(0.77,0.35,0.92,0.4); //with  writing about SF
leg2->SetFillColor(0);
leg2->SetBorderSize(0);
leg2->SetTextFont(42);
leg2->SetTextSize(0.025);

const int nhistos = 1; //40//52
TString hist_names[nhistos]={"h_PtZ_wt_"};
std::array<int,100> LOGY_array = {};

TString hist_names_sum[nhistos]={};
TString sum_histos_names[nhistos]={};
std::string hist_names_sort[nhistos];
for (int i=0;i<nhistos;i++){
	hist_names_sort[i] = hist_names[i];
	hist_names_sum[i] = hist_names[i];
	hist_names_sum[i].Append("_sum");
	sum_histos_names[i] = hist_names[i];
	sum_histos_names[i].Append("_sum0");
}
Float_t discriminators[nhistos];
TString stacks_names[nhistos];
for (int i=0;i<nhistos;i++){
	stacks_names[i] = hist_names[i];
	stacks_names[i].Prepend("s");
}
TString output_names[nhistos];
for (int i=0;i<nhistos;i++){
	output_names[i] = hist_names[i];
	output_names[i].Prepend(dir_name);
//	output_names[i].Prepend("triggercorr/");
//	output_names[i].Prepend("Aftertriggercorr2/");
	output_names[i].Append(set[set_type]);
//	output_names[i].Append("_amc_BDTcut.pdf");
	output_names[i].Append("_amc.pdf");
//	output_names[i].Append(set[set_type]+"_v25.png");
}

TH1F *data_histos[nhistos];
TH1F *data_histos2[nhistos];
TH1F *data_histosQCDup[nhistos];
TH1F *data_histosQCDlo[nhistos];
TH1F *data_histosJESup[nhistos];
TH1F *data_histosJESlo[nhistos];
TH1F *data_histosTrig[nhistos];
TH1F *signal_histos[nhistos];
TH1F *signal_histos2[nhistos];
TH1F *tthbb_histos[nhistos];
TH1F *tthnbb_histos[nhistos];
TH1F *gf_histos[nhistos];
TH1F *sum_histos[nhistos];
TH1F *sum_histosUp[nhistos];
TH1F *sum_histosDown[nhistos];
TH1F *histos_forUnc[nhistos];
TH1F *histos_for_legened[nhistos];
TH1F *discr_histos[nhistos];//Mqq,delta eta, delta phi, qgl, btag //12,13,14,21,22
TH1F *hBkgVis[nhistos];
TH1F *hBkgUncUp[nhistos];
TH1F *hBkgUncLo[nhistos];
TH1F *hBkgUncUpTrig[nhistos];
TH1F *hBkgUncLoTrig[nhistos];
TH1F *hBkgQCDUp[nhistos];
TH1F *hBkgQCDLo[nhistos];
TH1F *hBkgJESUp[nhistos];
TH1F *hBkgJESLo[nhistos];


int files=0; 
THStack *stacks[nhistos];
for (int i=0;i<nhistos;++i){
	stacks[i] = new THStack(stacks_names[i],"");
}
Double_t totalBG=0.;
Double_t totalQCD=0.;
Double_t totalMC=0.;
Double_t totalData=0.;
Double_t totalDataQCD=0.;
ofstream out_efficiency;
ofstream out_discrimination;
out_efficiency.open(dir_name+"efficiency.txt"); 




float ratio[5];
/*
int counter_ratio=0;
TString file_names_mc[nfiles];
do{
float dataInt = 0;
float MCint = 0;
float DYint =0;
file_names_mc[0] = file_names[0];
for (int i=1;i<nfiles;i++) {
	if (counter_ratio==0) file_names_mc[i] = file_names[i];
	if (counter_ratio==1) file_names_mc[i] = file_names_QCDup[i];
	if (counter_ratio==2) file_names_mc[i] = file_names_QCDdown[i];
	if (counter_ratio==3) file_names_mc[i] = file_names_JESup[i];
	if (counter_ratio==4) file_names_mc[i] = file_names_JESdown[i];
	if (((counter_ratio==1) || (counter_ratio==2))&& (i==5)) file_names_mc[i]=file_names[i];
}
files=0;
do{
	cout<<"here2"<<endl;
	TFile *file_initial_mc;		
	cout<<file_names_mc[files]<<endl;
  	file_initial_mc = TFile::Open(file_names_mc[files]);
	string file_name_tag = file_names_mc[files].Data();
	TH1F *histos_mc[100];
	int hist=0;
	histos_mc[hist] = (TH1F*)file_initial_mc->Get(hist_names[hist])->Clone("mc");
	if (files!=0 )histos_mc[hist]->Scale(lumi);  
	if (files==0) dataInt = histos_mc[hist]->Integral(0,histos_mc[hist]->GetNbinsX()+1); 
	if (files!=0) {
		if (file_name_tag.find("DYJetstoLL")!=std::string::npos)  DYint= histos_mc[hist]->Integral(0,histos_mc[hist]->GetNbinsX()+1);
		else MCint+=histos_mc[hist]->Integral(0,histos_mc[hist]->GetNbinsX()+1); 
	}
	files++;	
} while (files<nfiles); 
float tmp_ratio = (dataInt-MCint)/DYint;
cout<<tmp_ratio <<endl;
ratio[counter_ratio] = tmp_ratio;
//cout<<tmp_ratio<<endl;
counter_ratio++;
} while (counter_ratio<5);

*/
files=0;

do{
	TFile *file_initial;
  	file_initial = TFile::Open(file_names[files]);
	string file_name_tag = file_names[files].Data();
	string leg_name_tag = leg_names[files].Data();
	TH1F *histos[100];
	for (int hist=0;hist<nhistos;++hist){
		TString tmp_hist_tmp = "h_PtZ_wt_";
		tmp_hist_tmp.Append(file_hist_names[files]);
		hist_names[hist] = tmp_hist_tmp;
		histos[hist] = (TH1F*)file_initial->Get(hist_names[hist])->Clone("h");
		if (files==0) data_histos[hist] = (TH1F*)file_initial->Get(hist_names[hist])->Clone("data");
	//	if (files>0)histos[hist]->Scale(lumi); 
	//	if (file_name_tag.find("DYJetstoLL")!=std::string::npos)  histos[hist]->Scale(ratio[0]);  
//		if (file_name_tag.find("DYJetstoLL")!=std::string::npos)  histos[hist]->Scale(0.992);  
//		if (file_name_tag.find("EWK_LLJJ")!=std::string::npos)	histos[hist]->Scale(0.97957); //new signal x-sec
//		if (file_name_tag.find("EWK_LLJJ")!=std::string::npos)	histos[hist]->Scale(1.017); //new signal x-sec

		
////////////////////////////
////////////////////////////
////////////////////////////
// Signal now is included in the sum total MC
		if (files==1) 	sum_histos[hist] = (TH1F*)histos[hist]->Clone(sum_histos_names[hist]);
		if (files>1)	sum_histos[hist]->Add(histos[hist]); 

////////////////////////////
////////////////////////////
////////////////////////////
		if (files>0) histos[hist]->Sumw2(kFALSE);
//		if (hist==1) cout<<files<<"   "<<histos[1]->Integral() <<endl;

	//	if (files==1) {
		if (files==nfiles-1) {
			signal_histos[hist] = (TH1F*)file_initial->Get(hist_names[hist])->Clone(hist_names_sum[hist]+"newhist");
	//		signal_histos[hist]->Scale(lumi);
	//		signal_histos[hist]->Scale(0.97957); //new signal x-sec
		//	signal_histos[hist]->Scale(10);
			signal_histos[hist]->Sumw2(kFALSE);
		//	signal_histos[hist]->SetLineColor(LINECOLOR[files]);
			signal_histos[hist]->SetLineColor(kRed+2);
			signal_histos[hist]->SetLineStyle(LINESTYLE[files]);
			signal_histos[hist]->SetLineWidth(3);
			signal_histos2[hist]=(TH1F*)signal_histos[hist]->Clone("signalHist2");
		}
		histos[hist]->SetLineColor(LINECOLOR[files]);
		histos[hist]->SetLineStyle(LINESTYLE[files]);
		histos[hist]->SetLineWidth(LINEWIDTH[files]);
		histos[hist]->SetFillStyle(FILLSTYLE[files]);
		if ((files!=0)) histos[hist]->SetFillColor(FILLCOLOR[files]);
		
		if (files==0) {
			histos[hist]->SetMarkerStyle(20);
			data_histos[hist]->SetLineColor(1);
			data_histos[hist]->SetMarkerStyle(10);
			data_histos[hist]->SetMarkerSize(.8);
			data_histos[hist]->SetLineWidth(2);
		}
	 	//if (files>=bg_begin) stacks[hist]->Add(histos[hist]);
	 	if (files>=1) stacks[hist]->Add(histos[hist]);
		if (hist==0) histos_for_legened[files] = (TH1F*)histos[0]->Clone("newd");
		if (files==bg_begin)	discr_histos[hist] = (TH1F*)file_initial->Get(hist_names[hist])->Clone("discr");
		if ((files>bg_begin)&&(files<(nfiles-1))){
			discr_histos[hist]->Add(histos[hist]); 
   	}
		}
		if (files>=bg_begin) totalBG+=histos[0]->Integral(0,histos[0]->GetNbinsX()+1);
		if (files>0) totalMC+=histos[0]->Integral(0,histos[0]->GetNbinsX()+1);
		if (files==0) {totalData+=histos[0]->Integral(0,histos[0]->GetNbinsX()+1); }
		if (files==0) out_efficiency<<"Sample  \t\t\t yield(per "<< lumi<<" pb^-1)"<<endl;
		if (files==0) out_efficiency<<leg_names[order[files]]<<"\t \t \t"<< std::setprecision(8)<<histos[0]->Integral(0,histos[0]->GetNbinsX()+1) <<endl;
		else out_efficiency<<leg_names[order[files]]<<"\t\t\t  "<<std::setprecision(8)<<histos[0]->Integral(0,histos[0]->GetNbinsX()+1)<<endl;
		if (files==nfiles-1) out_efficiency<<"Total BG"<<"\t \t \t  "<<std::setprecision(8)<<totalBG<<endl;
		if (files==nfiles-1) out_efficiency<<"Total MC"<<"\t \t \t  "<<std::setprecision(8)<<totalMC<<endl;
		if (files==nfiles-1) out_efficiency<<"Data/MC"<<"\t \t \t  "<<std::setprecision(3)<<totalData/totalMC<<endl;
	//	if (files==nfiles-1) out_efficiency<<"DY MC scaled with "<<"\t \t \t  "<<std::setprecision(3)<<ratio[0]<<endl;
files++;
}while (files<nfiles);
out_efficiency.close();


/////////////////////////////////////
//aTGC files

TH1F *hist_aTGC[3];
for (int atgc = 0;atgc<3;atgc++){
	TString file_name_atgc;
//	if (atgc==0) {file_name_atgc = "new_files/histos_aTGC_plot_ptZ";
	if (atgc==0) {file_name_atgc = "histos_aTGC_plot_ptZ";
		file_name_atgc.Append(set[set_type]);
		file_name_atgc.Append("_SMfromATGC_cw0_cx6_cb750.root");
	}	
	if (atgc==1) {
	//	file_name_atgc = "new_files/histos_aTGC_plot_ptZ";
		file_name_atgc = "histos_aTGC_plot_ptZ";
		file_name_atgc.Append(set[set_type]);
		file_name_atgc.Append("_SMfromATGC_cwww0_cw20_cb750.root");
	}
//	if (atgc==2) {		file_name_atgc= "new_files/histos_aTGC_plot_ptZ";
	if (atgc==2) {		file_name_atgc= "histos_aTGC_plot_ptZ";
		file_name_atgc.Append(set[set_type]);
		file_name_atgc.Append("_SMfromATGC_cw0_cx6_cb750.root");
	}
	TString hist_name_atgc;
//	int tgc_color[3] = {kViolet+1,kMagenta,kBlue};
	int tgc_color[3] = {kBlue,kMagenta,kSpring};
	if (atgc==0) hist_name_atgc="aTGC_plot_cx_6_cb_0";
	if (atgc==1) hist_name_atgc="aTGC_plot_cw_20_cb_0";
	if (atgc==2) hist_name_atgc="aTGC_plot_cx_0_cb_750";
	TFile *file_atgc = TFile::Open(file_name_atgc);
	hist_aTGC[atgc] = (TH1F*)file_atgc->Get(hist_name_atgc);
	hist_aTGC[atgc]->SetLineWidth(3);
	hist_aTGC[atgc]->SetLineStyle(2+atgc);
	hist_aTGC[atgc]->SetLineColor(tgc_color[atgc]);
}






////////////////
files=3;
do{
	TFile *file_initial_up;		
   file_initial_up = TFile::Open(file_names[files]);
	TH1F *histos_QCDup[100];
	for (int hist=0;hist<nhistos;++hist){
		TString tmp_hist_tmp = "h_PtZ_wt_";
		tmp_hist_tmp.Append(file_hist_names[files]);
		tmp_hist_tmp.Append("_LHEscaleUp");
		hist_names[hist] = tmp_hist_tmp;
	//	cout<<hist_names[hist]<<endl;
		histos_QCDup[hist] = (TH1F*)file_initial_up->Get(hist_names[hist])->Clone("hup");
//		cout<<histos_QCDup[hist]->Integral()<<endl;
		histos_QCDup[hist]->SetLineWidth(2);
		if (files==3) 	hBkgQCDUp[hist] = (TH1F*)histos_QCDup[hist]->Clone(sum_histos_names[hist]+"QCDup");
		if (files>3)	hBkgQCDUp[hist]->Add(histos_QCDup[hist]);
		hBkgQCDUp[hist]->SetLineColor(kRed);
		hBkgQCDUp[hist]->SetLineStyle(2);
	}
	TFile *file_initial_down;
  	file_initial_down = TFile::Open(file_names[files]);
	TH1F *histos_QCDdown[100];
	for (int hist=0;hist<nhistos;++hist){
		TString tmp_hist_tmp = "h_PtZ_wt_";
		tmp_hist_tmp.Append(file_hist_names[files]);
		tmp_hist_tmp.Append("_LHEscaleDown");
		hist_names[hist] = tmp_hist_tmp;
		histos_QCDdown[hist] = (TH1F*)file_initial_down->Get(hist_names[hist])->Clone("hdown");
		histos_QCDdown[hist]->SetLineWidth(2);
		if (files==3) 	hBkgQCDLo[hist] = (TH1F*)histos_QCDdown[hist]->Clone(sum_histos_names[hist]+"QCDdown");
		if (files>3)	hBkgQCDLo[hist]->Add(histos_QCDdown[hist]);
		hBkgQCDLo[hist]->SetLineColor(kBlue);
		hBkgQCDLo[hist]->SetLineStyle(3);
	}
	files++;
}while (files<nfiles);

files=3;
do{
	TFile *file_initial_up;		
  	file_initial_up = TFile::Open(file_names[files]);
	TH1F *histos_JESup[100];
	for (int hist=0;hist<nhistos;++hist){
		TString tmp_hist_tmp = "h_PtZ_wt_";
		tmp_hist_tmp.Append(file_hist_names[files]);
		tmp_hist_tmp.Append("_JESUp");
		hist_names[hist] = tmp_hist_tmp;
		histos_JESup[hist] = (TH1F*)file_initial_up->Get(hist_names[hist])->Clone("hup");
		histos_JESup[hist]->SetLineWidth(2);
		if (files==3) 	hBkgJESUp[hist] = (TH1F*)histos_JESup[hist]->Clone(sum_histos_names[hist]+"JESup");
		if (files>3)	hBkgJESUp[hist]->Add(histos_JESup[hist]);
		hBkgJESUp[hist]->SetLineColor(kRed);
		hBkgJESUp[hist]->SetLineStyle(2);
	}
	TFile *file_initial_down;
  	file_initial_down = TFile::Open(file_names[files]);
	TH1F *histos_JESdown[100];
	for (int hist=0;hist<nhistos;++hist){
		TString tmp_hist_tmp = "h_PtZ_wt_";
		tmp_hist_tmp.Append(file_hist_names[files]);
		tmp_hist_tmp.Append("_JESDown");
		hist_names[hist] = tmp_hist_tmp;
		histos_JESdown[hist] = (TH1F*)file_initial_down->Get(hist_names[hist])->Clone("hdown");
		histos_JESdown[hist]->SetLineWidth(2);
		if (files==3) 	hBkgJESLo[hist] = (TH1F*)histos_JESdown[hist]->Clone(sum_histos_names[hist]+"JESdown");
		if (files>3)	hBkgJESLo[hist]->Add(histos_JESdown[hist]);
		hBkgJESLo[hist]->SetLineColor(kBlue);
		hBkgJESLo[hist]->SetLineStyle(3);
	}
	files++;
}while (files<nfiles);

/////////////////////////////////////


for (int hist=0;hist<nhistos;hist++){
	hBkgUncUp[hist] = (TH1F*)sum_histos[hist]->Clone("hBkgUncUp");
	hBkgUncLo[hist] = (TH1F*)sum_histos[hist]->Clone("hBkgUncLo");
  	hBkgVis[hist]   = (TH1F*)sum_histos[hist]->Clone("hbkgVis");
  	for(int i=0;i<hBkgUncUp[hist]->GetNbinsX();i++) {
  		float e = 0.0;
    	if (sum_histos[hist]->GetBinContent(i+1) != 0) {
      	e = sum_histos[hist]->GetBinError(i+1)/sum_histos[hist]->GetBinContent(i+1);
   	}
    	hBkgUncUp[hist]->SetBinContent(i+1,e);
    	hBkgUncLo[hist]->SetBinContent(i+1,-e);
	}
  	hBkgVis[hist]->SetMarkerSize(0);
 	hBkgVis[hist]->SetFillColor(kBlack);
 	hBkgVis[hist]->SetFillStyle(3004);
 	hBkgUncUp[hist]->SetLineColor(kBlack);
 	hBkgUncUp[hist]->SetLineWidth(1);
 	hBkgUncUp[hist]->SetFillColor(kBlack);
 	hBkgUncUp[hist]->SetFillStyle(3004);
	hBkgUncLo[hist]->SetLineColor(kBlack);
 	hBkgUncLo[hist]->SetLineWidth(1);
  	hBkgUncLo[hist]->SetFillColor(kBlack);
	hBkgUncLo[hist]->SetFillStyle(3004);
}

Float_t TSF[2] = {1.,1.};
Float_t kfactors[2] = {0.,0.};
//kfactors[set_type] = 1./(MC_data[set_type] * TSF[set_type]);
//kfactors[set_type] = 1./(ratio[0] * TSF[set_type]);
//cout<<"kfactor = "<<kfactors[set_type]<<endl;
//cout<<"Data/MC = "<< MC_data[set_type]<<endl;


histos_for_legened[0]->SetMarkerStyle(10);
histos_for_legened[0]->SetMarkerColor(1);
histos_for_legened[0]->SetLineColor(1);
histos_for_legened[0]->SetMarkerSize(0.8);
leg->AddEntry(histos_for_legened[0],leg_names[0],"LPE");
	//if (i==1)  leg->AddEntry(histos_for_legened[order_legend[i]],leg_names[i],"L");
//	if (i>=bg_begin) leg->AddEntry(histos_for_legened[order_legend[i]],leg_names[i],"F");

leg->AddEntry(histos_for_legened[1],leg_names_pas[1],"F");
leg->AddEntry(histos_for_legened[2],leg_names_pas[2],"F");
leg->AddEntry(histos_for_legened[3],leg_names_pas[3],"F");
leg->AddEntry(histos_for_legened[4],leg_names_pas[4],"F");
leg->AddEntry(signal_histos[0],"EW Zjj","L");
//leg->AddEntry(hist_aTGC[0],"aTGC #it{c_{www}} = 6","L");
//leg->AddEntry(hist_aTGC[1],"aTGC #it{c_{w}} = 20","L");
leg->AddEntry(hist_aTGC[0],"c_{#it{www}}/#Lambda^{2} = 6 TeV^{-2}","L");
leg->AddEntry(hist_aTGC[1],"c_{#it{w}}/#Lambda^{2} = 20 TeV^{-2}","L");
//leg->AddEntry(hist_aTGC[2],"aTGC c_{b} = 750","L");
leg->AddEntry(hBkgUncUp[0],"MC stat. unc.","F");



for (int d=0;d<nhistos;d++){
	discriminators[d] = Discr(discr_histos[d],signal_histos2[d]);
}

bubblesort(discriminators, hist_names_sort,nhistos);

//out_discrimination.open("Aftertriggercorr2/"+trigger[set_type]+dir_name+"discrimination.txt");
out_discrimination.open(dir_name+"discrimination.txt");
for (int d=0;d<nhistos;d++){
	if (d==0) out_discrimination<<"Variable &\t d"<<endl;
	out_discrimination<<"$"<<hist_names_sort[d]<<"$"<<" & \t "<< std::setprecision(2)<< discriminators[d]<<endl;
}
out_discrimination.close();

	float left2 = gStyle->GetPadLeftMargin();
	float right2 = gStyle->GetPadRightMargin();
	float top2 = gStyle->GetPadTopMargin();
	float bottom2 = gStyle->GetPadBottomMargin();
//TLatex* tex = new TLatex(0.75,0.95,"35.9 fb^{-1} (13 TeV)");
TLatex* tex = new TLatex(0.95,.95,"35.9 fb^{-1} (13 TeV)");
tex->SetNDC();
tex->SetTextAlign(35);
tex->SetTextFont(42);
tex->SetTextSize(0.035);
tex->SetLineWidth(2);
//TLatex *tex1 = new TLatex(0.17,0.95,"CMS");
TLatex *tex1 = new TLatex(0.21,0.95,"CMS");
tex1->SetNDC();
tex1->SetTextAlign(20);
tex1->SetTextFont(61);
tex1->SetTextSize(0.04);
tex1->SetLineWidth(2);
//TLatex* tex2 = new TLatex(0.27,0.89,"Work in progress");
TLatex* tex2 = new TLatex(0.26,0.885,"Preliminary");
tex2->SetNDC();
tex2->SetTextAlign(20);
tex2->SetTextFont(52);
tex2->SetTextSize(0.035);
tex2->SetLineWidth(2);
TString temp_str;
temp_str.Form("%2.2f",kfactors[set_type]);
//		temp_str.Form("%2.2f",(1./MC_data[set_type]));
TString k_factor_str;
TLatex* tex_set = new TLatex(0.667,0.86,set_names[set_type]);
//TLatex* tex_set = new TLatex(0.65,0.86,set_names[set_type]); ///for linear qgl
tex_set->SetNDC();
tex_set->SetTextAlign(20);
tex_set->SetTextFont(42);
tex_set->SetTextSize(0.03);
tex_set->SetLineWidth(2);
TLatex* tex_k = new TLatex(0.63,0.89,k_factor_str);
tex_k->SetNDC();
tex_k->SetTextAlign(20);
tex_k->SetTextFont(42);
tex_k->SetTextSize(0.03);
tex_k->SetLineWidth(2);
//	TPaveText pCMSset(0.57,1.-top2*2.,0.67,0.92,"NDC");
	TPaveText pCMSset(0.72,1.-top2*2.,0.9,0.92,"NDC");
	pCMSset.SetTextFont(42);
	pCMSset.SetTextSize(top2*0.75);
	pCMSset.SetTextAlign(12);
	pCMSset.SetFillStyle(-1);
	pCMSset.SetBorderSize(0);
	pCMSset.AddText(set_names[set_type]);
	TPaveText pBDT(0.4,1.-top2*2.,0.6,0.92,"NDC");
	pBDT.SetTextFont(42);
	pBDT.SetTextSize(top2*0.75);
	pBDT.SetTextAlign(12);
	pBDT.SetFillStyle(-1);
	pBDT.SetBorderSize(0);
	pBDT.AddText("BDT > 0.92");



	

for (int i=0;i<nhistos;i++){
		
		TString can_name="c1";
		TCanvas *c1 = new TCanvas(can_name,"",900,900);
		c1->cd();
		gPad->SetLogy(0);
		c1->SetBottomMargin(.3);
		bool LOGY=true;

	
	//	Double_t xmin = signal_histos[i]->GetBinCenter(0);
		Double_t xmin = 0;
	//	Double_t xmax =( signal_histos[i]->GetBinCenter((signal_histos[i]->GetNbinsX())-1) + signal_histos[i]->GetBinWidth((signal_histos[i]->GetNbinsX())-1)/2. );
		//Double_t xmax = signal_histos[i]->GetBinCenter((signal_histos[i]->GetNbinsX()))+signal_histos[i]->GetBinWidth((signal_histos[i]->GetNbinsX()));
		Double_t xmax = data_histos[i]->GetBinCenter((data_histos[i]->GetNbinsX()))+data_histos[i]->GetBinWidth((data_histos[i]->GetNbinsX()))/2.;
		cout<<xmin<<"  "<<xmax<<"  "<<data_histos[i]->GetBinError(data_histos[i]->GetNbinsX())<<endl;
//		data_histos[i]->SetBinError(data_histos[i]->GetNbinsX() ,TMath::Sqrt(data_histos[i]->GetBinContent(data_histos[i]->GetNbinsX())));
		string tmp_hist_name = hist_names[i].Data();
		TH1F *frame = new TH1F("frame","",1,xmin,xmax);
		TGaxis::SetExponentOffset(-0.07,0,"xy");
		frame->Reset();
		frame->SetMinimum(0.5);
      frame->SetMaximum(std::max(data_histos[i]->GetMaximum()*1.2,sum_histos[i]->GetMaximum()*1.2) );
		if (LOGY==true) {
			gPad->SetLogy();	
      //	frame->SetMaximum(std::max(data_histos[i]->GetMaximum()*200,sum_histos[i]->GetMaximum()*200) );
      	frame->SetMaximum(std::max(data_histos[i]->GetMaximum()*10,sum_histos[i]->GetMaximum()*10) );//for BDT cut
      	if (set_type==1)frame->SetMaximum(std::max(data_histos[i]->GetMaximum()*40,sum_histos[i]->GetMaximum()*40) );//for BDT cut
		}
		TGaxis::SetMaxDigits(4);
      frame->GetXaxis()->SetTitleOffset(0.91);
      frame->SetStats(0);
		frame->GetYaxis()->SetNdivisions(505);
	 	frame->GetXaxis()->SetLabelSize(0.0);
		char name[1000];
		string tmp_xaxis_name = signal_histos[i]->GetXaxis()->GetTitle();
      	sprintf(name,"Events / %1.0f %s",data_histos[i]->GetBinWidth(1),"GeV");
		frame->GetYaxis()->SetTitle(name);

      frame->Draw();
		tex->Draw();
		tex1->Draw();
	//	tex2->Draw();
	pCMSset.Draw("same");
//		if (tmp_hist_name.find("_bdt")!=std::string::npos)  pBDT.Draw("same");
	//	pBDT.Draw("same");
		leg->Draw("same");
    	stacks[i]->Draw("same");	
		signal_histos[i]->Draw("same");
		for (int atgc=0;atgc<2;atgc++)
			hist_aTGC[atgc]->Draw("same")	;	

/////////////////cross check of stupid THStack//////
//
//		sum_histos[i]->SetLineColor(kCyan);
//		sum_histos[i]->Scale(lumi);
//		sum_histos[i]->Draw("Lsame");
//
///////////////////////////////////////////////////
		data_histos[i]->Draw("PEsame");
		hBkgVis[i]->Draw("same E2");

	//	leg2->Draw("same");
  		gPad->RedrawAxis();
//	pCMSset.Draw("same");
	
  		TPad* pad2 = new TPad("pad2", "pad2", 0., 0., 1., 1.);
 		pad2->SetTopMargin(0.73);
 	//	pad2->SetRightMargin(0.05);  //0.25
 		pad2->SetFillColor(0);
  		pad2->SetFillStyle(0);
  		pad2->Draw();
  		pad2->cd(0);
  //		gPad->SetGridy();

		TH1F *frame2 = new TH1F("frame2","",1,xmin,xmax);
	//	frame2->SetMinimum(-.5);
    //  frame2->SetMaximum(.5);
		frame2->SetMinimum(-1.);//for bdt
      frame2->SetMaximum(1.);//fot bdt
	//	if ((hist_names[i].CompareTo("hAdJetHT_bdt")==0)||(hist_names[i].CompareTo("hJet3_pt_bdt")==0)) {
		if (tmp_hist_name.find("_bdt")!=std::string::npos) {
			frame2->SetMinimum(-1.);
    	 	frame2->SetMaximum(1.);
		}
      frame2->SetStats(0);
      frame2->SetTitleFont(42,"x");
		frame2->SetTitleFont(42,"y");
      frame2->SetTitleSize(0.13, "XYZ");
		frame2->GetYaxis()->SetNdivisions(505);
 		frame2->GetYaxis()->SetTickLength(0.06);
  		frame2->GetYaxis()->SetTitleSize(0.04);
  		frame2->GetYaxis()->SetTitleOffset(1.5);
 		frame2->GetYaxis()->SetLabelSize(0.03);
  		frame2->GetYaxis()->CenterTitle(kTRUE);
  		frame2->GetXaxis()->SetTitleSize(0.05);
  		frame2->GetXaxis()->SetLabelSize(0.04);
		frame2->SetXTitle("p_{TZ} (GeV)");
		frame2->SetYTitle("Data / MC - 1");
		frame2->Draw();	

 		Double_t aa[2] = {xmin,xmax};
   	Double_t bb[2] = {0,0};
   	TGraph *cons = new TGraph(2,aa,bb);
    	cons->SetLineStyle(2);
		cons->Draw("Lsame");
		

		data_histos2[i] = (TH1F*)data_histos[i]->Clone("new");
		data_histosQCDup[i] = (TH1F*)hBkgQCDUp[i]->Clone("newQCDup");
		data_histosQCDup[i]->SetLineColor(kBlue);
		data_histosQCDup[i]->SetLineStyle(2);
		data_histosQCDlo[i] = (TH1F*)hBkgQCDLo[i]->Clone("newQCDlow");
		data_histosQCDlo[i]->SetLineStyle(2);
		data_histosQCDlo[i]->SetLineColor(kBlue);
		data_histosQCDlo[i]->Scale(data_histos[i]->Integral()/data_histosQCDlo[i]->Integral());
		data_histosQCDup[i]->Scale(data_histos[i]->Integral()/data_histosQCDup[i]->Integral());
		cout<<data_histosQCDlo[i]->Integral()<< "    "<<data_histosQCDup[i]->Integral()<<endl;
		data_histosJESup[i] = (TH1F*)hBkgJESUp[i]->Clone("newJESup");
		data_histosJESup[i]->SetLineColor(kRed);
		data_histosJESup[i]->SetLineStyle(2);
		data_histosJESlo[i] = (TH1F*)hBkgJESLo[i]->Clone("newJESlow");
		data_histosJESlo[i]->SetLineStyle(2);
		data_histosJESlo[i]->SetLineColor(kRed);
		data_histosJESlo[i]->Scale(data_histos[i]->Integral()/data_histosJESlo[i]->Integral());
		data_histosJESup[i]->Scale(data_histos[i]->Integral()/data_histosJESup[i]->Integral());
		cout<<data_histosJESlo[i]->Integral()<< "    "<<data_histosJESup[i]->Integral()<<endl;
	//	float chi2_data_mc = data_histos2[i]->Chi2Test(sum_histos[i],"UWCHI2/NDF");
		float ks = data_histos2[i]->KolmogorovTest(sum_histos[i]);
		data_histos2[i]->Add(sum_histos[i],-1);
		data_histos2[i]->Divide(sum_histos[i]);
		for (int j=0;j<data_histos2[i]->GetNbinsX();j++){
			if (sum_histos[i]->GetBinContent(j+1)!=0) data_histos2[i]->SetBinError(j+1, TMath::Sqrt(data_histos[i]->GetBinContent(j+1))/sum_histos[i]->GetBinContent(j+1)) ;
			else if (data_histos[i]->GetBinContent(j+1)!=0) data_histos2[i]->SetBinError(j+1, TMath::Sqrt(data_histos[i]->GetBinContent(j+1))/data_histos[i]->GetBinContent(j+1)) ;
					else data_histos2[i]->SetBinError(j+1,0.); 
		}
		data_histos2[i]->Draw("PEsame");
		data_histosQCDup[i]->Add(sum_histos[i],-1);
		data_histosQCDup[i]->Divide(sum_histos[i]);
		data_histosQCDup[i]->Draw("HISTsame");
		data_histosQCDlo[i]->Add(sum_histos[i],-1);
		data_histosQCDlo[i]->Divide(sum_histos[i]);
		data_histosQCDlo[i]->Draw("HISTsame");
		//cout<<data_histosQCDlo[i]->Integral()<< "    "<<data_histosQCDup[i]->Integral()<<endl;
		data_histosJESup[i]->Add(sum_histos[i],-1);
		data_histosJESup[i]->Divide(sum_histos[i]);
		data_histosJESup[i]->Draw("HISTsame");
		data_histosJESlo[i]->Add(sum_histos[i],-1);
		data_histosJESlo[i]->Divide(sum_histos[i]);
		data_histosJESlo[i]->Draw("HISTsame");
		hBkgUncUp[i]->Draw("HIST same");
		hBkgUncLo[i]->Draw("HIST same");
		
//		temp_str.Form("%2.2f",chi2_data_mc);
//		TString chi2_str = temp_str.Prepend("#chi^{2} = ");
//		TLatex *chi2_latex = new TLatex(0.21,0.23,temp_str);
 //	  	chi2_latex->SetNDC();
  	// 	chi2_latex->SetTextAlign(20);
 // 		chi2_latex->SetTextFont(42);
  // 	chi2_latex->SetTextSize(0.025);
  // 	chi2_latex->SetLineWidth(2);
	//	chi2_latex->Draw();
		temp_str.Form("%2.2f",ks);
		TString ks_str = temp_str.Prepend(", KS = ");
		TLatex *ks_latex = new TLatex(0.31,0.23,temp_str);
 	  	ks_latex->SetNDC();
  	 	ks_latex->SetTextAlign(20);
  		ks_latex->SetTextFont(42);
   	ks_latex->SetTextSize(0.025);
   	ks_latex->SetLineWidth(2);
	//	ks_latex->Draw();
		TLegend *leg_jes = new TLegend(0.16,0.13,0.3,0.29); //without writing about SF
//		TLegend *leg_jes = new TLegend(0.16,0.13,0.3,0.26); //without writing about SF
		leg_jes->SetFillStyle(0);
		leg_jes->SetBorderSize(0);
		leg_jes->SetTextFont(42);
		leg_jes->SetTextSize(0.025);
	//	leg_jes->AddEntry(data_histosQCDup[i],"QCD scale up/down","L");
		leg_jes->AddEntry(data_histosQCDup[i],"#mu_{F}, #mu_{R} scale up/down","L");
		leg_jes->AddEntry(data_histosJESup[i],"JES up/down","L");
		leg_jes->Draw();
		
		cout<<output_names[i]<<endl;

		pad2->RedrawAxis();
		c1->Print(output_names[i]);
		c1->Delete();
	}



return 0;
}
