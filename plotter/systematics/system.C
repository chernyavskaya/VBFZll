#include <ctype.h>
#include <iomanip>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <cstdlib>
#include "TFile.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2.h"
#include "TTree.h"
#include "TChain.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <algorithm> 
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
#include <TProfile.h>
#include "math.h"
#include "/afs/cern.ch/work/n/nchernya/VBFZll/plotter/EWcorr.C"
#include "/mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/mucorr/2016/RoccoR.cc"
//#include "EWcorr.C"

Double_t erf( Double_t *x, Double_t *par){
  return par[0]/2.*(1.+TMath::Erf((x[0]-par[1])/par[2]));
}
Double_t erf2( Double_t *x, Double_t *par){
  return par[0]/2.*(1.+TMath::Erf((x[0]-par[1])/par[2]))+ (1.-par[0]);
}


Double_t inter_func(Double_t *x, Double_t *par){
	Double_t a0 = par[0];
		Double_t a1 = par[1];
		Double_t a2 = par[2];
		Double_t a3 = par[3];
		Double_t a4 = par[4];
		Double_t a5 = par[5];
		Double_t a6 = par[6];
		Double_t a7 = par[7];
		Double_t y = a0 + a1*x[0] + a2*pow(x[0],2)+ a3*pow(x[0],3)+ a4*pow(x[0],4)+ a5*pow(x[0],5)+ a6*pow(x[0],6)+ a7*pow(x[0],7);
		if (y<0)  y=0;
	return y; 
}
Double_t inter_ratio(Double_t *x, Double_t *par){
		Double_t a0 = par[0];
		Double_t a1 = par[1];
		Double_t a2 = par[2];
		Double_t a3 = par[3];
		Double_t a4 = par[4];
		Double_t a5 = par[5];
		Double_t a6 = par[6];
		Double_t a7 = par[7];
		Double_t a02 = par[8];
		Double_t a12 = par[9];
		Double_t a22 = par[10];
		Double_t a32 = par[11];
		Double_t a42 = par[12];
		Double_t a52 = par[13];
		Double_t a62 = par[14];
		Double_t a72 = par[15];
		Double_t y1 = a0 + a1*x[0] + a2*pow(x[0],2)+ a3*pow(x[0],3)+ a4*pow(x[0],4)+ a5*pow(x[0],5)+ a6*pow(x[0],6)+ a7*pow(x[0],7);
		Double_t y2 = a02 + a12*x[0] + a22*pow(x[0],2)+ a32*pow(x[0],3)+ a42*pow(x[0],4)+ a52*pow(x[0],5)+ a62*pow(x[0],6)+ a72*pow(x[0],7);
		Double_t y=0;
	if (y2<0) y2=0;
	if (y1<0)  y1=0;
	if (y2!=0) 	y = y1*y1/y2;
	else  y=0;
	return y;
}
		







#define SWAP2(A, B) { TLorentzVector t = A; A = B; B = t; }
void SortByEta(std::vector<TLorentzVector> &jets){
  int i, j;
	int n=jets.size();
  for (i = n - 1; i >= 0; i--){
    for (j = 0; j < i; j++){
      if (jets[j].Eta() < jets[j + 1].Eta() ){
        SWAP2( jets[j], jets[j + 1] );
		}
    }
	}
}

float getScaleFactor(TH2F *scaleMap, double pt, double eta, float sf_err, bool abs) {
//    std::cout<<"called getScaleFactor"<<std::endl;
  //  std::cout<<pt<<":, "<<eta<<std::endl;
    float sfactor = 1.0;
    int binx = scaleMap->GetXaxis()->FindBin(pt);
	 int biny;
    if (abs==0) biny = scaleMap->GetYaxis()->FindBin(eta);
    else biny = scaleMap->GetYaxis()->FindBin(TMath::Abs(eta));
  //  std::cout<<binx<<": ,"<<biny<<std::endl;
    if ( (binx != 0) && (binx != scaleMap->GetNbinsX()+1) && (biny != 0) && (biny != scaleMap->GetNbinsY()+1)) {
        sfactor = scaleMap->GetBinContent(binx, biny);
        sf_err = scaleMap->GetBinError(binx, biny);
	//		cout<<sfactor<<endl;
        if (sfactor == 0.0) {
            // bin was not filled for w/e reason, assume we don't have value in this 2D bin from the json
            sfactor = 1.0;
            sf_err = 0.0;
        }
    }
    //std::cout<<sfactor<<std::endl;
    return sfactor;
}


float getScaleFactor1D(TH1F *scaleMap, double eta, float sf_err, bool abs) {
//    std::cout<<"called getScaleFactor"<<std::endl;
  //  std::cout<<pt<<":, "<<eta<<std::endl;
    float sfactor = 1.0;
	 int binx;
    if (abs==0) binx = scaleMap->GetXaxis()->FindBin(eta);
    else binx = scaleMap->GetXaxis()->FindBin(TMath::Abs(eta));
    if ( (binx != 0) && (binx != scaleMap->GetNbinsX()+1) ) {
        sfactor = scaleMap->GetBinContent(binx);
        sf_err = scaleMap->GetBinError(binx);
	//		cout<<sfactor<<endl;
        if (sfactor == 0.0) {
            // bin was not filled for w/e reason, assume we don't have value in this 2D bin from the json
            sfactor = 1.0;
            sf_err = 0.0;
        }
    }
    //std::cout<<sfactor<<std::endl;
    return sfactor;
}



typedef std::map<double, int> JetList;
const int njets = 30;




typedef struct {
   Float_t eta[njets];
   Float_t pt[njets];
   Float_t JEC_corr[njets];
   Float_t JEC_corr_up[njets];
   Float_t JEC_corr_down[njets];
   Float_t JER_corr[njets];
   Float_t JER_corr_up[njets];
   Float_t JER_corr_down[njets];
   Float_t phi[njets];
	Float_t mass[njets];
	Float_t btag[njets];
	Int_t nsoft;
	Float_t soft_pt[njets];
	Float_t soft_eta[njets];
	Float_t soft_mass[njets];
	Float_t qgl[njets];
	Int_t nsoft2;
	Int_t nsoft5;
	Int_t nsoft10;
	Int_t EWKnsoft2;
	Int_t EWKnsoft5;
	Int_t EWKnsoft10;
	Int_t id[njets];
	Int_t puId[njets];
	Float_t HTsoft;
	Int_t partonFlavour[njets];
	Float_t EWKHTsoft;
	Float_t pt_regVBF[njets];	
	Float_t ptd[njets];
	Float_t axis2[njets];
	Int_t mult[njets];
	Float_t leadTrackPt[njets];
	Float_t blike_VBF[njets];
} Jets;

using namespace std;


int main(int argc, char* argv[]){



TString file_name = std::string(argv[1]);
TString file_tag = std::string(argv[2]);
TString region = std::string(argv[3]); 
int data = atoi(argv[4]);
TString heppyVersion = std::string(argv[5]);
TString postfix = std::string(argv[6]);
TString output = std::string(argv[7]);

std::map <TString, float> xsec;
xsec["SingleMuon"] = 1.;
xsec["SingleMuonB"] = 1.;
xsec["SingleMuonC"] = 1.;
xsec["SingleMuonD"] = 1.;
xsec["SingleMuonE"] = 1.;
xsec["SingleMuonF"] = 1.;
xsec["SingleMuonG"] = 1.;
xsec["SingleElectron"] =  1.;
xsec["SingleElectronB"] =  1.;
xsec["SingleElectronC"] =  1.;
xsec["SingleElectronD"] =  1.;
xsec["SingleElectronE"] =  1.;
xsec["SingleElectronF"] =  1.;
xsec["SingleElectronG"] =  1.;
xsec["DYJetstoLL"] =  5765.4;
xsec["DYJetstoLL_amc"] =  5765.4;
xsec["DYJetstoLL_HT100"] =  5765.4;
//xsec["DYJetstoLL_HT100_200"] = 173.96106;
//xsec["DYJetstoLL_HT200_400"] = 48.27802 ;
//xsec["DYJetstoLL_HT400_600"] =6.68755 ;
//xsec["DYJetstoLL_HT600_Inf"] = 2.588804;
//xsec["DYJetstoLL_HT100_200"] = 147.40 ; 
//xsec["DYJetstoLL_HT200_400"] = 40.99 ; 
//xsec["DYJetstoLL_HT400_600"] = 5.678 ; 
//xsec["DYJetstoLL_HT600_Inf"] = 2.198; 
xsec["DYJetstoLL_HT100_200"] = 181.302; 
xsec["DYJetstoLL_HT200_400"] =50.4177  ; 
xsec["DYJetstoLL_HT400_600"] =6.98394; 
xsec["DYJetstoLL_HT600_Inf"] =2.70354 ;
xsec["DYJetstoLL_HT600_800"] = 1.6814;
xsec["DYJetstoLL_HT800_1200"] = 0.7754;
xsec["DYJetstoLL_HT1200_2500"] = 0.186;
xsec["DYJetstoLL_HT2500_Inf"] = 0.00438495;
 
xsec["DYJetstoLL_Pt-100_amc"] = 5765.4; 
xsec["DYJetstoLL_Pt-100To250_amc"] = 83.12; 
xsec["DYJetstoLL_Pt-250To400_amc"] =3.047 ; 
xsec["DYJetstoLL_Pt-400To650_amc"] = 0.3921 ; 
xsec["DYJetstoLL_Pt-650ToInf_amc"] = 0.03636 ;
 
xsec["DYJetstoLL_amc_0J"] = 4585.27; //4732.;  normalize to 5765.4 pb, 1.032 
xsec["DYJetstoLL_amc_1J"] = 853.198;//880.5; 
xsec["DYJetstoLL_amc_2J"] = 325.194;//335.6; 
xsec["DYJetstoLL_amc_2J_ext"] = 325.194;//335.6; 


xsec["TT"] =809.;
xsec["WW"] =118.7;
xsec["WZ"] = 47.13;
xsec["ZZ"] =16.523;
xsec["EWK_LLJJ"]=1.630;
xsec["EWK_LLJJ_herwig"]=1.630;
xsec["interference"]=1.630;

xsec["QCD_HT100to200"] = 27990000;
xsec["QCD_HT200to300"] = 1712000 ;
xsec["QCD_HT300to500"] = 347700;
xsec["QCD_HT500to700"] = 32100 ;
xsec["QCD_HT700to1000"] = 6831;
xsec["QCD_HT1000to1500"] = 1207 ;
xsec["QCD_HT1500to2000"] =  119.9;
xsec["QCD_HT2000toInf"] = 25.24;

xsec["ST_tW"] = 71.7 ;			//inclusive decays
xsec["ST_tW_top"] = 35.85  ;			//inclusive decays
xsec["ST_tW_antitop"] = 35.85  ;			//inclusive decays
xsec["ST_s-channel"] = 3.36; //leptonic decays
xsec["ST_t-channel_top_4f_inclusiveDecays"] = 136.02;
xsec["ST_t-channel_antitop_4f_inclusiveDecays"] = 80.95;
xsec["ST_t-channel_top_4f_leptonDecays"] = 44.33;  //leptonDecays  , multiplied with BR 0.325
xsec["ST_t-channel_antitop_4f_leptonDecays"] = 26.38;//leptonDecays ,multiplied with BR 0.325

xsec["WJetsToLNu_amc"]  = 61526.7; //not going to use these ones
xsec["WJetsToLNu"]  = 61526.7;
xsec["TTZToLLNuNu"] = 0.2529;
xsec["tZq_ll"]=0.0758;

//k factors 1.21 are not included// 
/*xsec["WJetsToLNu_HT100To200"] = 1345 ;
xsec["WJetsToLNu_HT200To400"] = 359.7  ;
xsec["WJetsToLNu_HT400To600"] = 48.91;
xsec["WJetsToLNu_HT600To800"] =12.05;
xsec["WJetsToLNu_HT800To1200"] = 5.501;
xsec["WJetsToLNu_HT1200To2500"] = 1.329;
xsec["WJetsToLNu_HT2500ToInf"] = 0.03216;
*/
xsec["WJetsToLNu_HT100"]  = 61526.7;
xsec["WJetsToLNu_HT100To200"] = 1627.45 ;
xsec["WJetsToLNu_HT200To400"] = 435.236  ;
xsec["WJetsToLNu_HT400To600"] = 59.18109;
xsec["WJetsToLNu_HT600To800"] =14.5805;
xsec["WJetsToLNu_HT800To1200"] = 6.656210;
xsec["WJetsToLNu_HT1200To2500"] = 1.608089;
xsec["WJetsToLNu_HT2500ToInf"] = 0.0389135;
xsec["TTZToLLNuNu"] = 0.2529;
xsec["tZq_ll"]=0.0758;

xsec["EWKinterference"]=0.1701;



 int counter=0;



    
float gen_pos=0; 
float gen_neg=0; 
float gen_pos_weight=0; 
float gen_neg_weight=0; 


	
	Float_t presel=0;
	Float_t presel_vtype[10] = {0,0,0,0,0,0,0,0,0};
	Float_t pos_puweight=0;
	Float_t all_puweight=0.;
	Float_t puweight;
	Float_t puweightUp;
	Float_t puweightDown;
	Float_t PU=1.;
	Float_t genweight;
	Float_t bTagWeight;
	Float_t genweight0;
	float  trigWeight_tree;
	Int_t global_counter = 0;
	Int_t HLT_QuadPFJet_DoubleBTag_CSV_VBF_Mqq200;
	Int_t HLT_QuadPFJet_SingleBTag_CSV_VBF_Mqq460;
	Int_t HLT_IsoMu22;
	Int_t HLT_IsoTkMu22;
	Int_t HLT_IsoMu24;
	Int_t HLT_IsoTkMu24;
	Int_t HLT_IsoMu27;
	Int_t HLT_IsoTkMu27;
	Int_t HLT_Ele27_eta2p1;
	TFile *file_initial;
	TChain *tree_initial;

	Int_t nvLeptons, nselLeptons;
	const int brLeptons=13;
	Float_t vLeptons_pt[30], vLeptons_eta[30], vLeptons_phi[30], vLeptons_mass[30], vLeptons_SF_IdCutLoose[30], vLeptons_SF_IdCutTight[30], vLeptons_SF_IsoLoose[30], vLeptons_SF_IsoTight[30],vLeptons_SF_trk_eta[30], vLeptons_SF_HLT_RunD4p2[30],vLeptons_SF_HLT_RunD4p3[30] ;
	Int_t vLeptons_charge[30], vLeptons_pdgId[30],vLeptons_trackerLayers[30]; 
	TString str_leptons[brLeptons] = {"vLeptons_pt", "vLeptons_eta", "vLeptons_phi", "vLeptons_mass", "vLeptons_charge", "vLeptons_pdgId", "vLeptons_SF_IdCutLoose", "vLeptons_SF_IdCutTight", "vLeptons_SF_IsoLoose","vLeptons_SF_IsoTight","vLeptons_SF_trk_eta","vLeptons_SF_HLT_RunD4p2","vLeptons_SF_HLT_RunD4p3"};
	Float_t selLeptons_pt[30], selLeptons_eta[30], selLeptons_phi[30], selLeptons_mass[30], selLeptons_SF_IdCutLoose[30], selLeptons_SF_IdCutTight[30], selLeptons_SF_IsoLoose[30], selLeptons_SF_IsoTight[30],selLeptons_SF_trk_eta[30], selLeptons_SF_HLT_RunD4p2[30],selLeptons_SF_HLT_RunD4p3[30], selLeptons_relIso04[30] ;
	Int_t selLeptons_charge[30], selLeptons_pdgId[30], selLeptons_looseIdPOG[30], selLeptons_trackerLayers[30],  selLeptons_eleMVAIdSppring16GenPurp[30]; 

//	float interference_kfactor = 12./17.;
	float interference_kfactor = 1.;


	TFile* file_id_mu_bf = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/v25/TriggerEffMap_MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta_RunBCDEF.root");
	TH2F* id_mu_bf = (TH2F*)file_id_mu_bf->Get("TriggerEffMap_MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta");
	TFile* file_id_mu_aft = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/v25/TriggerEffMap_MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta_RunGH.root");
	TH2F* id_mu_aft = (TH2F*)file_id_mu_aft->Get("TriggerEffMap_MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta");

	TFile* file_trig_mu_bf = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/v25/TriggerEffMap_IsoMu24_OR_IsoTkMu24_PtEtaBins_RunBCDEF.root");
	TH2F* trig_mu_bf = (TH2F*)file_trig_mu_bf->Get("TriggerEffMap_IsoMu24_OR_IsoTkMu24_PtEtaBins");
	TFile* file_trig_mu_aft = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/v25/TriggerEffMap_IsoMu24_OR_IsoTkMu24_PtEtaBins_RunGH.root");
	TH2F* trig_mu_aft = (TH2F*)file_trig_mu_aft->Get("TriggerEffMap_IsoMu24_OR_IsoTkMu24_PtEtaBins");

	TFile* file_iso_mu_bf = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/v25/TriggerEffMap_LooseISO_LooseID_pt_eta_RunBCDEF.root");
	TH2F* iso_mu_bf = (TH2F*)file_iso_mu_bf->Get("TriggerEffMap_LooseISO_LooseID_pt_eta");
	TFile* file_iso_mu_aft = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/v25/TriggerEffMap_LooseISO_LooseID_pt_eta_RunGH.root");
	TH2F* iso_mu_aft = (TH2F*)file_iso_mu_aft->Get("TriggerEffMap_LooseISO_LooseID_pt_eta");

	TFile* file_track_mu_bf = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/v25/TriggerEffMap_Muons_trk_SF_RunBCDEF.root");
	TH1F* track_mu_bf = (TH1F*)file_track_mu_bf->Get("TriggerEffMap_Graph");
	TFile* file_track_mu_aft = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/v25/TriggerEffMap_Muons_trk_SF_RunGH.root");
	TH1F* track_mu_aft = (TH1F*)file_track_mu_aft->Get("TriggerEffMap_Graph");

	TFile* file_tracker_el = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/v25/TriggerEffMap_ScaleFactor_tracker_80x.root");
	TH2F* tracker_el = (TH2F*)file_tracker_el->Get("TriggerEffMap_ScaleFactor_tracker_80x");
	TFile* file_trig_el = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/v25/TriggerEffMap_Tight27AfterIDISO.root");
	TH2F* trig_el = (TH2F*)file_trig_el->Get("TriggerEffMap_Tight27AfterIDISO");
	TFile* file_id_el = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/v25/TriggerEffMap_EIDISO_ZH.root");
	TH2F* id_el = (TH2F*)file_id_el->Get("TriggerEffMap_EIDISO_ZH");



	
////////////////////////////

	
	file_initial = TFile::Open(file_name);
	
	tree_initial = (TChain*)file_initial->Get("tree");
	Float_t events_generated;
	Float_t events_generated_muF_QCDUp;
	Float_t events_generated_muF_QCDDown;
	Float_t events_generated_muR_QCDUp;
	Float_t events_generated_muR_QCDDown;
	Float_t events_generated_muFR_QCDUp;
	Float_t events_generated_muFR_QCDDown;
	TH1F *countPos;
	TH1F *countNeg;
	TH1F *countLHEScale;
	TH1F *countLHEPdf;
	TH1F *countWeighted;
	if ((data!=1)){
		countPos = (TH1F*)file_initial->Get("CountPosWeight");
		countNeg = (TH1F*)file_initial->Get("CountNegWeight");
 		countWeighted = (TH1F*)file_initial->Get("CountWeighted");
 		countLHEScale = (TH1F*)file_initial->Get("CountWeightedLHEWeightScale");
		countLHEPdf=	(TH1F*)file_initial->Get("CountWeightedLHEWeightPdf");
 		events_generated = countWeighted->GetBinContent(1);
		if (events_generated==0) events_generated =  countPos->GetBinContent(1) - countNeg->GetBinContent(1);
		events_generated_muF_QCDUp = countLHEScale->GetBinContent( countLHEScale->FindBin(0) );
		events_generated_muF_QCDDown = countLHEScale->GetBinContent( countLHEScale->FindBin(1) );
		events_generated_muR_QCDUp = countLHEScale->GetBinContent( countLHEScale->FindBin(2) );
		events_generated_muR_QCDDown = countLHEScale->GetBinContent( countLHEScale->FindBin(3) );
		events_generated_muFR_QCDUp = countLHEScale->GetBinContent( countLHEScale->FindBin(4) );
		events_generated_muFR_QCDDown = countLHEScale->GetBinContent( countLHEScale->FindBin(5) );
	} else events_generated = 1;
/*	if (file_tag.CompareTo("EWK_LLJJ")==0)  {
		events_generated=events_generated/2.00703;
		events_generated_muF_QCDUp=events_generated/2.00703;
		events_generated_muF_QCDDown=events_generated/2.00703;
		events_generated_muR_QCDUp=events_generated/2.00703;
		events_generated_muR_QCDDown=events_generated/2.00703;
		events_generated_muFR_QCDUp=events_generated/2.00703;
		events_generated_muFR_QCDDown=events_generated/2.00703;

	}*/
    Jets Jet;
    Float_t v_type;
    Float_t wrong_type=0.;
    Int_t nJets;
	Float_t JSON;	
	Float_t nPVs;	
	Float_t rho;	
	Float_t lheHT;	
	Float_t lheV_pt;	

	const int Nsyst = 10;
	Float_t bdt;
	int passSel;
	int passSel_JESR[4];
	Float_t bdt_JESR[4];
	Float_t met_pt;
	Float_t met_phi;

	Jets GenHiggsSisters;

	int pos_weight_presel=0;
 	Int_t selLeptons_tightId[20];
	Float_t selLeptons_relIso03[20] , selLeptons_chargedHadRelIso03[20], selLeptons_pfRelIso03[20];
	Float_t vLeptons_dz[20], vLeptons_edz[20];

	Int_t nGenVbosons;
	Float_t GenVbosons_pt[1];
	Int_t GenVbosons_pdgId[1];
	Float_t VtypeSim; 

Float_t LHE_weights_pdf_wgt[103];
Float_t LHE_weights_pdf_eigen[100];
Float_t LHE_weights_scale_wgt[10];
	
	float V_mass;
	ULong64_t evt;

    tree_initial->SetBranchAddress("Vtype",&v_type);
    tree_initial->SetBranchAddress("V_mass",&V_mass);
    tree_initial->SetBranchAddress("rho",&rho);
    tree_initial->SetBranchAddress("nJet",&nJets);
    tree_initial->SetBranchAddress("Jet_pt",Jet.pt);
    tree_initial->SetBranchAddress("Jet_corr_JECUp",Jet.JEC_corr_up);
    tree_initial->SetBranchAddress("Jet_corr_JECDown",Jet.JEC_corr_down);
    tree_initial->SetBranchAddress("Jet_corr",Jet.JEC_corr);
    tree_initial->SetBranchAddress("Jet_corr_JERUp",Jet.JER_corr_up);
    tree_initial->SetBranchAddress("Jet_corr_JERDown",Jet.JER_corr_down);
    tree_initial->SetBranchAddress("Jet_corr_JER",Jet.JER_corr);
    tree_initial->SetBranchAddress("Jet_eta",Jet.eta);
    tree_initial->SetBranchAddress("Jet_phi",Jet.phi);
	tree_initial->SetBranchAddress("Jet_mass",Jet.mass);
	tree_initial->SetBranchAddress("Jet_btagCSV",Jet.btag);
	tree_initial->SetBranchAddress("Jet_blike_VBF",Jet.blike_VBF);
	tree_initial->SetBranchAddress("Jet_id",Jet.id);	
	tree_initial->SetBranchAddress("Jet_puId",Jet.puId);
 	tree_initial->SetBranchAddress("Jet_leadTrackPt",Jet.leadTrackPt);
 	tree_initial->SetBranchAddress("Jet_partonFlavour",Jet.partonFlavour);
	tree_initial->SetBranchAddress("met_pt",&met_pt);
	tree_initial->SetBranchAddress("met_phi",&met_phi);
	
	tree_initial->SetBranchAddress("softActivityJets_pt",Jet.soft_pt);
	tree_initial->SetBranchAddress("softActivityJets_eta",Jet.soft_eta);
	tree_initial->SetBranchAddress("softActivityJets_mass",Jet.soft_mass);
	tree_initial->SetBranchAddress("softActivity_HT",&Jet.HTsoft);
	tree_initial->SetBranchAddress("softActivity_njets2",&Jet.nsoft2);
	tree_initial->SetBranchAddress("softActivity_njets5",&Jet.nsoft5);
	tree_initial->SetBranchAddress("softActivity_njets10",&Jet.nsoft10);
	tree_initial->SetBranchAddress("softActivityEWK_HT",&Jet.EWKHTsoft);
	tree_initial->SetBranchAddress("softActivityEWK_njets2",&Jet.EWKnsoft2);
	tree_initial->SetBranchAddress("softActivityEWK_njets5",&Jet.EWKnsoft5);
	tree_initial->SetBranchAddress("softActivityEWK_njets10",&Jet.EWKnsoft10);
	tree_initial->SetBranchAddress("Jet_qgl",Jet.qgl);
	tree_initial->SetBranchAddress("genWeight",&genweight);
	tree_initial->SetBranchAddress("bTagWeight",&bTagWeight);
	tree_initial->SetBranchAddress("puWeight",&puweight);
	tree_initial->SetBranchAddress("puWeightUp",&puweightUp);
	tree_initial->SetBranchAddress("puWeightDown",&puweightDown);
	tree_initial->SetBranchAddress("nPVs",&nPVs);
	tree_initial->SetBranchAddress("Jet_ptd",Jet.ptd);
	tree_initial->SetBranchAddress("Jet_axis2",Jet.axis2);
	tree_initial->SetBranchAddress("Jet_mult",Jet.mult);
 	tree_initial->SetBranchAddress("HLT_BIT_HLT_QuadPFJet_DoubleBTagCSV_VBF_Mqq200_v",&HLT_QuadPFJet_DoubleBTag_CSV_VBF_Mqq200);
 	tree_initial->SetBranchAddress("HLT_BIT_HLT_QuadPFJet_SingleBTagCSV_VBF_Mqq460_v",&HLT_QuadPFJet_SingleBTag_CSV_VBF_Mqq460);
 	tree_initial->SetBranchAddress("HLT_BIT_HLT_IsoMu22_v",&HLT_IsoMu22);
 	tree_initial->SetBranchAddress("HLT_BIT_HLT_IsoTkMu22_v",&HLT_IsoTkMu22);
 	tree_initial->SetBranchAddress("HLT_BIT_HLT_IsoMu24_v",&HLT_IsoMu24);
 	tree_initial->SetBranchAddress("HLT_BIT_HLT_IsoTkMu24_v",&HLT_IsoTkMu24);
 	tree_initial->SetBranchAddress("HLT_BIT_HLT_IsoMu27_v",&HLT_IsoMu27);
 	tree_initial->SetBranchAddress("HLT_BIT_HLT_IsoMu24_v",&HLT_IsoMu24);
 	tree_initial->SetBranchAddress("HLT_BIT_HLT_IsoTkMu27_v",&HLT_IsoTkMu27);
 	tree_initial->SetBranchAddress("HLT_BIT_HLT_Ele27_eta2p1_WPTight_Gsf_v",&HLT_Ele27_eta2p1);
	tree_initial->SetBranchAddress("Jet_pt_regVBF",Jet.pt_regVBF);
	tree_initial->SetBranchAddress("json",&JSON);
    
	tree_initial->SetBranchAddress("GenHiggsSisters_pt",GenHiggsSisters.pt);
    tree_initial->SetBranchAddress("GenHiggsSisters_eta",GenHiggsSisters.eta);
    tree_initial->SetBranchAddress("GenHiggsSisters_phi",GenHiggsSisters.phi);
	tree_initial->SetBranchAddress("GenHiggsSisters_mass",GenHiggsSisters.mass);
	tree_initial->SetBranchAddress("nGenVbosons",&nGenVbosons);
	tree_initial->SetBranchAddress("GenVbosons_pt",GenVbosons_pt);
	tree_initial->SetBranchAddress("GenVbosons_pdgId",GenVbosons_pdgId);
	tree_initial->SetBranchAddress("VtypeSim",&VtypeSim);

	tree_initial->SetBranchAddress("nvLeptons",&nvLeptons);
	tree_initial->SetBranchAddress("vLeptons_pt",vLeptons_pt);
	tree_initial->SetBranchAddress("vLeptons_eta",vLeptons_eta);
	tree_initial->SetBranchAddress("vLeptons_phi",vLeptons_phi);
	tree_initial->SetBranchAddress("vLeptons_mass",vLeptons_mass);
	tree_initial->SetBranchAddress("vLeptons_charge",vLeptons_charge);
	tree_initial->SetBranchAddress("vLeptons_pdgId",vLeptons_pdgId);
	tree_initial->SetBranchAddress("vLeptons_SF_IdCutLoose",vLeptons_SF_IdCutLoose);
	tree_initial->SetBranchAddress("vLeptons_SF_IdCutTight",vLeptons_SF_IdCutTight);
	tree_initial->SetBranchAddress("vLeptons_SF_IsoLoose",vLeptons_SF_IsoLoose);
	tree_initial->SetBranchAddress("vLeptons_SF_IsoTight",vLeptons_SF_IsoTight);
	tree_initial->SetBranchAddress("vLeptons_SF_trk_eta",vLeptons_SF_trk_eta);
	tree_initial->SetBranchAddress("vLeptons_SF_HLT_RunD4p2",vLeptons_SF_HLT_RunD4p2);
	tree_initial->SetBranchAddress("vLeptons_SF_HLT_RunD4p3",vLeptons_SF_HLT_RunD4p3);
	tree_initial->SetBranchAddress("vLeptons_trackerLayers", vLeptons_trackerLayers);
	tree_initial->SetBranchAddress("nselLeptons",&nselLeptons);
	tree_initial->SetBranchAddress("selLeptons_pt",selLeptons_pt);
	tree_initial->SetBranchAddress("selLeptons_eta",selLeptons_eta);
	tree_initial->SetBranchAddress("selLeptons_phi",selLeptons_phi);
	tree_initial->SetBranchAddress("selLeptons_mass",selLeptons_mass);
	tree_initial->SetBranchAddress("selLeptons_charge",selLeptons_charge);
	tree_initial->SetBranchAddress("selLeptons_pdgId",selLeptons_pdgId);
	tree_initial->SetBranchAddress("selLeptons_looseIdPOG",selLeptons_looseIdPOG);
	tree_initial->SetBranchAddress("selLeptons_relIso04",selLeptons_relIso04);
	tree_initial->SetBranchAddress("selLeptons_relIso03",selLeptons_relIso03);
   tree_initial->SetBranchAddress("selLeptons_eleMVAIdSppring16GenPurp",selLeptons_eleMVAIdSppring16GenPurp); 
	tree_initial->SetBranchAddress("selLeptons_trackerLayers", selLeptons_trackerLayers);
	


//	for (int i=0;i<brLeptons;i++){
//		tree_initial->SetBranchAddress(str_leptons[i],);
//	}	

	tree_initial->SetBranchAddress("LHE_weights_pdf_wgt",LHE_weights_pdf_wgt);
	tree_initial->SetBranchAddress("LHE_weights_pdf_eigen",LHE_weights_pdf_eigen);
	tree_initial->SetBranchAddress("LHE_weights_scale_wgt",LHE_weights_scale_wgt);
	tree_initial->SetBranchAddress("lheHT",&lheHT);
	tree_initial->SetBranchAddress("lheV_pt",&lheV_pt);
	tree_initial->SetBranchAddress("BDT_VBF",&bdt);
	tree_initial->SetBranchAddress("BDT_VBF_JES_up",&bdt_JESR[0]);
	tree_initial->SetBranchAddress("BDT_VBF_JES_down",&bdt_JESR[1]);
	tree_initial->SetBranchAddress("BDT_VBF_JER_up",&bdt_JESR[2]);
	tree_initial->SetBranchAddress("BDT_VBF_JER_down",&bdt_JESR[3]);

	tree_initial->SetBranchAddress("PassSelection_JES_up",&passSel_JESR[0]);
	tree_initial->SetBranchAddress("PassSelection_JES_down",&passSel_JESR[1]);
	tree_initial->SetBranchAddress("PassSelection_JER_up",&passSel_JESR[2]);
	tree_initial->SetBranchAddress("PassSelection_JER_down",&passSel_JESR[3]);
	tree_initial->SetBranchAddress("PassSelection_nom",&passSel);


	tree_initial->SetBranchAddress("evt",&evt);

	if (data==1){
		genweight = 1.;
		bTagWeight = 1.;
		puweight=1.;
	}

	float lumi = 35900;


	TString uncertainty_name[Nsyst] = {"","puWeight","LHE_weights_scale","LHE_weights_scale_muF","LHE_weights_scale_muR","JES","JER","MDG_NLO_corr","int_shape","QGL"};
	TH1F *hbdtUp[20];
	TH1F *hbdt_atanhUp[20];
	TH1F *hbdtDown[20];
	TH1F *hbdt_atanhDown[20];
	TH1F *hbdt_atanhUp_pdf[20];
	TH1F *hbdt_atanhDown_pdf[20];
	for (int i=0;i<Nsyst;i++){
		char histTitleUp[150];
		char histTitleDown[150];
		if (i==0) {
			sprintf(histTitleUp,"BDT_\%s_\%s",region.Data(),file_tag.Data());
			hbdtUp[i] = new TH1F(histTitleUp,"",500,-1.,1.);
			sprintf(histTitleUp,"atanhBDT_\%s_\%s",region.Data(),file_tag.Data());
			hbdt_atanhUp[i] = new TH1F(histTitleUp,"",900,0.,3.);
		}
		else {
			sprintf(histTitleUp,"BDT_\%s_\%s_CMS_ewkzjj_\%sUp",region.Data(),file_tag.Data(),uncertainty_name[i].Data());
			sprintf(histTitleDown,"BDT_\%s_\%s_CMS_ewkzjj_\%sDown",region.Data(),file_tag.Data(),uncertainty_name[i].Data());
			hbdtUp[i] = new TH1F(histTitleUp,"",500,-1.,1.);
			hbdtDown[i] = new TH1F(histTitleDown,"",500,-1.,1.);
			sprintf(histTitleUp,"atanhBDT_\%s_\%s_CMS_ewkzjj_\%sUp",region.Data(),file_tag.Data(),uncertainty_name[i].Data());
			sprintf(histTitleDown,"atanhBDT_\%s_\%s_CMS_ewkzjj_\%sDown",region.Data(),file_tag.Data(),uncertainty_name[i].Data());
			hbdt_atanhUp[i] = new TH1F(histTitleUp,"",900,0.,3.);
			hbdt_atanhDown[i] = new TH1F(histTitleDown,"",900,0.,3.);
		}
		}	

	for (int i=0;i<Nsyst;i++){
		hbdt_atanhUp[i]->Sumw2();
		if (i!=0) {
			hbdt_atanhDown[i]->Sumw2();
			hbdtDown[i]->Sumw2();
		}
		hbdtUp[i]->Sumw2();
	}


	/*
	for (int i=0;i<30;i++){
		char histTitleUp[150];
		char histTitleDown[150];
		sprintf(histTitleUp,"atanhBDT_\%s_\%s_CMS_ewkzjj_pdf\%iUp",region.Data(),file_tag.Data(),i);
		sprintf(histTitleDown,"atanhBDT_\%s_\%s_CMS_ewkzjj_pdf\%iDown",region.Data(),file_tag.Data(),i);
		hbdt_atanhUp_pdf[i] = new TH1F(histTitleUp,"",900,0.,3.);
		hbdt_atanhDown_pdf[i] = new TH1F(histTitleDown,"",900,0.,3.);
	}
	*/

	int nentries = tree_initial->GetEntries() ;
	


	TF1* func_Mqq = new TF1("func_Mqq","pol6",0,10);
	TF1* func_Mqq_up = new TF1("func_Mqq_up","pol6",0,10);
	TF1* func_Mqq_down = new TF1("func_Mqq_down","pol6",0,10);
	func_Mqq->FixParameter(0,-1378.361758);
	func_Mqq->FixParameter(1,1327.383178);
	func_Mqq->FixParameter(2,-529.261759);
	func_Mqq->FixParameter(3,111.854413);
	func_Mqq->FixParameter(4,-13.208941);
	func_Mqq->FixParameter(5,0.826140);
	func_Mqq->FixParameter(6,-0.021376);
	

	func_Mqq_up->FixParameter(0,-935.385276);
	func_Mqq_down->FixParameter(0,-1821.330360);
	func_Mqq_up->FixParameter(1,904.196013);
	func_Mqq_down->FixParameter(1,1750.563094);
	func_Mqq_up->FixParameter(2,-361.222999);
	func_Mqq_down->FixParameter(2,-697.297753);
	func_Mqq_up->FixParameter(3,76.352839);
	func_Mqq_down->FixParameter(3,147.355426);
	func_Mqq_up->FixParameter(4,-8.999610);
	func_Mqq_down->FixParameter(4,-17.418208);
	func_Mqq_up->FixParameter(5,0.560532);
	func_Mqq_down->FixParameter(5,1.091744);
	func_Mqq_up->FixParameter(6,-0.014406);
	func_Mqq_down->FixParameter(6,-0.028346);

	TF1* interference_func = new TF1("inter_fnew",inter_func,0,10,8);
interference_func->FixParameter(0,-3236.73471);
interference_func->FixParameter(1,3158.70903);
interference_func->FixParameter(2,-1314.92713);
interference_func->FixParameter(3,302.849052);
interference_func->FixParameter(4,-41.69131);
interference_func->FixParameter(5,3.43119691);
interference_func->FixParameter(6,-0.15633696);
interference_func->FixParameter(7,0.00304252817);

//fixed scale
	TF1* interference_func_fixed = new TF1("inter_fold",inter_func,5,10,8);
interference_func_fixed->FixParameter(0,-2829.38504);
interference_func_fixed->FixParameter(1,2723.93909);
interference_func_fixed->FixParameter(2,-1118.42065);
interference_func_fixed->FixParameter(3,254.03378);
interference_func_fixed->FixParameter(4,-34.4852553);
interference_func_fixed->FixParameter(5,2.79853847);
interference_func_fixed->FixParameter(6,-0.125730001);
interference_func_fixed->FixParameter(7,0.00241279903);
	TF1* interference_down = new TF1("interference_down",inter_ratio,5,10,16);
interference_down->FixParameter(0,-3236.73471);
interference_down->FixParameter(1,3158.70903);
interference_down->FixParameter(2,-1314.92713);
interference_down->FixParameter(3,302.849052);
interference_down->FixParameter(4,-41.69131);
interference_down->FixParameter(5,3.43119691);
interference_down->FixParameter(6,-0.15633696);
interference_down->FixParameter(7,0.00304252817);
interference_down->FixParameter(8,-2829.38504);
interference_down->FixParameter(9,2723.93909);
interference_down->FixParameter(10,-1118.42065);
interference_down->FixParameter(11,254.03378);
interference_down->FixParameter(12,-34.4852553);
interference_down->FixParameter(13,2.79853847);
interference_down->FixParameter(14,-0.125730001);
interference_down->FixParameter(15,0.00241279903);



	TF1* func_qgl_q = new TF1("func_qgl_q","pol3",0.,1.);
	func_qgl_q->FixParameter(0,0.981581);
	func_qgl_q->FixParameter(1,-0.255505);
	func_qgl_q->FixParameter(2,0.929524);
	func_qgl_q->FixParameter(3,-0.666978);
	TF1* func_qgl_g = new TF1("func_qgl_g","pol7",0.,1.);
	func_qgl_g->FixParameter(0,0.612992);
	func_qgl_g->FixParameter(1,6.27);
	func_qgl_g->FixParameter(2,-34.3663);
	func_qgl_g->FixParameter(3,92.8668);
	func_qgl_g->FixParameter(4,-99.927);
	func_qgl_g->FixParameter(5,-21.1421);
	func_qgl_g->FixParameter(6, 113.218);
	func_qgl_g->FixParameter(7,-55.7067);
//pythia8 quark (|eta|<2.0, pT inclusive, pythia ): -0.666978*x*x*x + 0.929524*x*x -0.255505*x + 0.981581
//pythia8 gluon (|eta|<2.0, pT inclusive, pythia ): -55.7067*x^7 + 113.218*x^6 -21.1421*x^5 -99.927*x^4 + 92.8668*x^3 -34.3663*x^2 + 6.27*x + 0.612992



	RoccoR  *rc = new RoccoR("/mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/mucorr/2016/rcdata.2016.v3/");


TFile file(output+"/"+file_tag+"_"+region+"_"+heppyVersion+"_"+postfix+".root","recreate");
float scaleWeightsUp[20];
float scaleWeightsDown[20];
string file_tag_str = file_tag.Data();
TString uncNotAppl[20] = {"","data","WW_WZ_ZZ_ST_tW","WW_WZ_ZZ_ST_tW","WW_WZ_ZZ_ST_tW","","","",""};
for (int current_syst=0;current_syst<Nsyst;current_syst++){
//for (int current_syst=0;current_syst<1;current_syst++){
//	if ((current_syst!=0) && (current_syst!=9)) continue;
///////////////////////////////////////////////
//to save time don't calculate muF and muR only LHE scale unc, only combined/////////
	if ((current_syst==3) || (current_syst==4)) continue;
///////////////////////////////////////////////
	if ((data==1)&&(current_syst!=0)) continue;
	if (( file_tag.CompareTo("EWKinterference")==0 )&&(current_syst!=0)) continue;
	string uncNotAppl_str = uncNotAppl[current_syst].Data();
	if (current_syst!=0) if (uncNotAppl_str.find(file_tag_str)!=std::string::npos) continue;
	if  ( (uncertainty_name[current_syst].CompareTo("MDG_NLO_corr")==0) &&(  !((file_tag.CompareTo("DYJetstoLL")==0)  || (file_tag_str.find("DYJetstoLL_HT")!=std::string::npos) ) ))  continue;
	if  ( (uncertainty_name[current_syst].CompareTo("int_shape")==0) &&(  !((file_tag.CompareTo("interference")==0) ) ))  continue;
	bdt = 0;
	passSel=0;
	genweight = 0 ;
	cout<<current_syst<<endl;
	for (int entry=0; entry<nentries;++entry){
//	for (int entry=0; entry<1000;++entry){
        tree_initial->GetEntry(entry);

	
		if  ((uncertainty_name[current_syst].CompareTo("JER")!=0) &&  (uncertainty_name[current_syst].CompareTo("JES")!=0)) {
			if (passSel!=1) continue;
		}
	
	//	if ((file_tag.CompareTo("EWK_LLJJ")==0) &&(evt%2!=0)) continue;
	
		scaleWeightsUp[0] = puweight/events_generated;
		scaleWeightsDown[0] = puweight/events_generated;

		scaleWeightsUp[1] = puweightUp/events_generated;
		scaleWeightsDown[1] = puweightDown/events_generated;
		scaleWeightsUp[2] = puweight/events_generated_muFR_QCDUp*LHE_weights_scale_wgt[4];
		scaleWeightsDown[2] = puweight/events_generated_muFR_QCDDown*LHE_weights_scale_wgt[5];
		scaleWeightsUp[3] = puweight/events_generated_muF_QCDUp*LHE_weights_scale_wgt[0];
		scaleWeightsDown[3] = puweight/events_generated_muF_QCDDown*LHE_weights_scale_wgt[1];
		scaleWeightsUp[4] = puweight/events_generated_muR_QCDUp*LHE_weights_scale_wgt[2];
		scaleWeightsDown[4] = puweight/events_generated_muR_QCDDown*LHE_weights_scale_wgt[3];

		scaleWeightsUp[5] = puweight/events_generated;
		scaleWeightsDown[5] = puweight/events_generated;
		scaleWeightsUp[6] = puweight/events_generated;
		scaleWeightsDown[6] = puweight/events_generated;

		scaleWeightsUp[7] = puweight/events_generated;
		scaleWeightsDown[7] = puweight/events_generated;
		scaleWeightsUp[8] = puweight/events_generated;
		scaleWeightsDown[8] = puweight/events_generated;
		scaleWeightsUp[9] = puweight/events_generated;
		scaleWeightsDown[9] = puweight/events_generated;

		if (JSON!=1) {
			continue;
		}


//		if (region.CompareTo("mu")==0) if (!(v_type==0)) continue;
//		if (region.CompareTo("el")==0) if (!(v_type==1)) continue;

			

		if (data==1) PU=1.;
		else PU=puweight;
		genweight=genweight/TMath::Abs(genweight)*xsec[file_tag]*lumi;

		
		double GenVbosons_pt_first = GenVbosons_pt[0];
		int GenVbosons_pdgId_first = GenVbosons_pdgId[0];


 
		if  ((file_tag.CompareTo("DYJetstoLL_HT100")==0)) if (lheHT>100) continue;  
		if  ((file_tag.CompareTo("DYJetstoLL_Pt-100_amc")==0)) if (lheV_pt>100) continue;  

		int pt_num1 = -1;
		int pt_num2 = -1;
		TLorentzVector Qjet1;
		TLorentzVector Qjet2;
		TLorentzVector qq;
		int good_jets = 0;
		vector<TLorentzVector> jets_pv;
		vector<int> jets_indices;
		///////////////////////
		//preselection/////
		//////////////////////

		for (int i=0;i<nJets;i++){
			TLorentzVector jet0;
			if (!((Jet.id[i]>2)&&(Jet.puId[i]>0))) continue;
			jet0.SetPtEtaPhiM(Jet.pt[i],Jet.eta[i],Jet.phi[i],Jet.mass[i]);
			jets_pv.push_back(jet0);
			jets_indices.push_back(i);
			good_jets++;
		}
		Qjet1 = jets_pv[0];
		Qjet2 = jets_pv[1];
		float jet3_pt = jets_pv[2].Pt();
		qq=Qjet1+Qjet2;
		Float_t Mqq = qq.M();
		Float_t qq_pt = qq.Pt();
		Float_t qqDeltaEta = TMath::Abs(Qjet1.Eta()-Qjet2.Eta());
		Float_t qqDeltaPhi = TMath::Abs(Qjet1.DeltaPhi(Qjet2));
//////////////////leptons////////////////
		TLorentzVector lepton1;
		TLorentzVector lepton2;
		TLorentzVector Zll;
		int idx_1stLepton = 0;
		int idx_2ndLepton = 0;
		int count_l=0;
		if (region.CompareTo("el")==0) {
		for (int i=0; i<nselLeptons;i++ ){
			if (!((selLeptons_eleMVAIdSppring16GenPurp[i]>=2)&& (selLeptons_relIso03[i]<0.15)&& (TMath::Abs(selLeptons_pdgId[i])==11))) continue;
		//	if ((count_l==1) && (selLeptons_charge[idx_1stLepton]*selLeptons_charge[i] > 0)) continue;
			if (count_l==1) {
				idx_2ndLepton=i;
				lepton2.SetPtEtaPhiM(selLeptons_pt[idx_2ndLepton], selLeptons_eta[idx_2ndLepton], selLeptons_phi[idx_2ndLepton], selLeptons_mass[idx_2ndLepton]);
				count_l++;
				break;
			}
			if (count_l==0) {
				idx_1stLepton=i;
				lepton1.SetPtEtaPhiM(selLeptons_pt[idx_1stLepton], selLeptons_eta[idx_1stLepton], selLeptons_phi[idx_1stLepton], selLeptons_mass[idx_1stLepton]);
				count_l++;
			}
		}
		}
		if (region.CompareTo("mu")==0) {
		count_l=0;
		idx_1stLepton = 0;
		idx_2ndLepton = 0;
		for (int i=0; i<nselLeptons;i++ ){
			if (!((selLeptons_looseIdPOG[i]>0) && (selLeptons_relIso04[i]<0.25) && (TMath::Abs(selLeptons_pdgId[i])==13 )) ) continue;
		//	if ((count_l==1) && (selLeptons_charge[idx_1stLepton]*selLeptons_charge[i] > 0)) continue;
			if (count_l==1) {
				idx_2ndLepton=i;
				lepton2.SetPtEtaPhiM(selLeptons_pt[idx_2ndLepton], selLeptons_eta[idx_2ndLepton], selLeptons_phi[idx_2ndLepton], selLeptons_mass[idx_2ndLepton]);
				count_l++;
				break;
			}
			if (count_l==0) {
				idx_1stLepton=i;
				lepton1.SetPtEtaPhiM(selLeptons_pt[idx_1stLepton], selLeptons_eta[idx_1stLepton], selLeptons_phi[idx_1stLepton], selLeptons_mass[idx_1stLepton]);
				count_l++;
			}
		}
		}
		if (count_l<2)  continue;
		if ((selLeptons_charge[idx_1stLepton]*selLeptons_charge[idx_2ndLepton]) >0) continue;
///////////////muon corrections 2016 calibration////////////////
		if (region.CompareTo("mu")==0) {
			double dataSF1, dataSF2 ;
			double mcSF1, mcSF2;
			if (data==1) {
				dataSF1 = rc->kScaleDT(selLeptons_charge[idx_1stLepton], lepton1.Pt(), lepton1.Eta(), lepton1.Phi(), 0, 0);
				dataSF2 = rc->kScaleDT(selLeptons_charge[idx_2ndLepton], lepton2.Pt(), lepton2.Eta(), lepton2.Phi(), 0, 0);
				lepton1.SetPtEtaPhiM(lepton1.Pt()*dataSF1,lepton1.Eta(), lepton1.Phi(), lepton1.M() );
				lepton2.SetPtEtaPhiM(lepton2.Pt()*dataSF2,lepton2.Eta(), lepton2.Phi(), lepton2.M() );
			}
			if (data!=1) {
				double u1 = gRandom->Rndm();
				double u2 = gRandom->Rndm();
				mcSF1 = rc->kScaleAndSmearMC(selLeptons_charge[idx_1stLepton], lepton1.Pt(), lepton1.Eta(), lepton1.Phi(),  selLeptons_trackerLayers[idx_1stLepton], u1, u2, 0, 0);
				mcSF2 = rc->kScaleAndSmearMC(selLeptons_charge[idx_2ndLepton], lepton2.Pt(), lepton2.Eta(), lepton2.Phi(),  selLeptons_trackerLayers[idx_2ndLepton], u1, u2, 0, 0);
				lepton1.SetPtEtaPhiM(lepton1.Pt()*mcSF1,lepton1.Eta(), lepton1.Phi(), lepton1.M() );
				lepton2.SetPtEtaPhiM(lepton2.Pt()*mcSF2,lepton2.Eta(), lepton2.Phi(), lepton2.M() );
			}
		}
		

////////////////////////////////////////////////////////////////
			if  (region.CompareTo("mu")==0) if (!((HLT_IsoMu24==1) || (HLT_IsoTkMu24==1)  )) continue; 
			if  (region.CompareTo("el")==0) if (!(HLT_Ele27_eta2p1 == 1)) continue;

		if (data!=1) {
			if (region.CompareTo("mu")==0) {
				float SF_mu_bf_err1 = 0.;
				float SF_mu_bf_err2 = 0.;
				float SF_mu_aft_err1 = 0.;
				float SF_mu_aft_err2 = 0.;
				bool abs=1;
				float eff1 =20.1/36.4*getScaleFactor(trig_mu_bf, lepton1.Pt(), lepton1.Eta(), SF_mu_bf_err1,abs ) + 16.3/36.4*getScaleFactor(trig_mu_aft, lepton1.Pt(), lepton1.Eta(), SF_mu_bf_err1,abs ) ;  	
	
				float eff1_id =20.1/36.4*getScaleFactor(id_mu_bf, lepton1.Pt(), lepton1.Eta(), SF_mu_bf_err1,abs ) + 16.3/36.4*getScaleFactor(id_mu_aft, lepton1.Pt(), lepton1.Eta(), SF_mu_bf_err1,abs ) ;  	
				float eff2_id =20.1/36.4*getScaleFactor(id_mu_bf, lepton2.Pt(), lepton2.Eta(), SF_mu_bf_err2,abs ) + 16.3/36.4*getScaleFactor(trig_mu_aft, lepton2.Pt(), lepton2.Eta(), SF_mu_bf_err2,abs  ) ;  
				float eff1_iso =20.1/36.4*getScaleFactor(iso_mu_bf, lepton1.Pt(), lepton1.Eta(), SF_mu_bf_err1,abs ) + 16.3/36.4*getScaleFactor(id_mu_aft, lepton1.Pt(), lepton1.Eta(), SF_mu_bf_err1,abs ) ;  	
				float eff2_iso =20.1/36.4*getScaleFactor(iso_mu_bf, lepton2.Pt(), lepton2.Eta(), SF_mu_bf_err2,abs ) + 16.3/36.4*getScaleFactor(trig_mu_aft, lepton2.Pt(), lepton2.Eta(), SF_mu_bf_err2,abs  ) ;  
				abs=0; 
				float eff1_tracker =20.1/36.4*getScaleFactor1D(track_mu_bf, lepton1.Eta(), SF_mu_bf_err1,abs ) + 16.3/36.4*getScaleFactor1D(track_mu_aft,  lepton1.Eta(), SF_mu_bf_err1,abs );  	
				float eff2_tracker =20.1/36.4*getScaleFactor1D(track_mu_bf, lepton2.Eta(), SF_mu_bf_err1,abs ) + 16.3/36.4*getScaleFactor1D(track_mu_aft,  lepton2.Eta(), SF_mu_bf_err1,abs );  	

				genweight*= eff1*eff1_id*eff2_id*eff1_iso*eff2_iso*eff1_tracker*eff2_tracker; 	
			}
			if (region.CompareTo("el")==0) {
				float SF_el_err1 = 0.;
				float SF_el_err2 = 0.;
				bool abs=0;
				float eff1 = getScaleFactor(trig_el, lepton1.Pt(), lepton1.Eta(), SF_el_err1,abs );  	
				float eff1_id =getScaleFactor(id_el, lepton1.Pt(), lepton1.Eta(), SF_el_err1,abs ) ;  	
				float eff2_id =getScaleFactor(id_el, lepton2.Pt(), lepton2.Eta(), SF_el_err2,abs  ) ;  
				float eff1_tracker =getScaleFactor(tracker_el, lepton1.Pt(), lepton1.Eta(), SF_el_err1,abs ) ;  	
				float eff2_tracker =getScaleFactor(tracker_el, lepton2.Pt(), lepton2.Eta(), SF_el_err2,abs  ) ;  
				genweight*=eff1* eff1_id*eff2_id*eff1_tracker*eff2_tracker; 	
			}
		}

		string file_tag_str = file_tag.Data();
		if  (file_tag_str.find("DYJetstoLL")!=std::string::npos)  genweight*=ptWeightEWK(nGenVbosons, GenVbosons_pt_first, VtypeSim, GenVbosons_pdgId_first); 

		float alpha_Mqq_up=0.;
		float alpha_Mqq_down=0.;
		
		Float_t Mqq_log = TMath::Log(Mqq);	
	//	if ((Mqq_log< 8.3227 )&&(Mqq_log>5.101)) if ( (file_tag_str.find("DYJetstoLL_HT")!=std::string::npos) || (file_tag.CompareTo("DYJetstoLL")==0)) {		
		if ( (file_tag_str.find("DYJetstoLL_HT")!=std::string::npos) || (file_tag.CompareTo("DYJetstoLL")==0)) {
				if (uncertainty_name[current_syst].CompareTo("MDG_NLO_corr")!=0) { 
					genweight*=func_Mqq->Eval(Mqq_log);
				} else if (uncertainty_name[current_syst].CompareTo("MDG_NLO_corr")==0){
					alpha_Mqq_up = func_Mqq_up->Eval(Mqq_log);
					alpha_Mqq_down = func_Mqq_down->Eval(Mqq_log);
				}
		}


		float qgl_weight=1.;
		if (!( (data==1)|| (Jet.partonFlavour[jets_indices[0]] ==0 ) || (TMath::Abs(Jet.eta[jets_indices[0]])>=2) || (Jet.qgl[jets_indices[0]] < 0) ) ) {
			if (TMath::Abs(Jet.partonFlavour[jets_indices[0]]) < 4 ) qgl_weight*=func_qgl_q->Eval(Jet.qgl[jets_indices[0]]);
			if (TMath::Abs(Jet.partonFlavour[jets_indices[0]]) ==21 ) qgl_weight*=func_qgl_g->Eval(Jet.qgl[jets_indices[0]]);
		}
	//	cout<<qgl_weight<<endl;
		if (!( (data==1)|| (Jet.partonFlavour[jets_indices[1]] ==0 ) || (TMath::Abs(Jet.eta[jets_indices[1]])>=2) || (Jet.qgl[jets_indices[1]] < 0)) ) {
			if (TMath::Abs(Jet.partonFlavour[jets_indices[1]]) < 4 ) qgl_weight*=func_qgl_q->Eval(Jet.qgl[jets_indices[1]]);
			if (TMath::Abs(Jet.partonFlavour[jets_indices[1]]) ==21 ) qgl_weight*=func_qgl_g->Eval(Jet.qgl[jets_indices[1]]);
		}
			
		scaleWeightsUp[0]*= qgl_weight;//nominal is with puweight
		scaleWeightsUp[1]*= qgl_weight;//nominal is with puweight
		scaleWeightsUp[2]*= qgl_weight;//nominal is with puweight
		scaleWeightsUp[3]*= qgl_weight;//nominal is with puweight
		scaleWeightsUp[4]*= qgl_weight;//nominal is with puweight
		scaleWeightsUp[5]*= qgl_weight;//nominal is with puweight
		scaleWeightsUp[6]*= qgl_weight;//nominal is with puweight
		scaleWeightsUp[7]*= qgl_weight;//nominal is with puweight
		scaleWeightsUp[8]*= qgl_weight;//nominal is with puweight
		scaleWeightsDown[0]*= qgl_weight;//nominal is with puweight
		scaleWeightsDown[1]*= qgl_weight;//nominal is with puweight
		scaleWeightsDown[2]*= qgl_weight;//nominal is with puweight
		scaleWeightsDown[3]*= qgl_weight;//nominal is with puweight
		scaleWeightsDown[4]*= qgl_weight;//nominal is with puweight
		scaleWeightsDown[5]*= qgl_weight;//nominal is with puweight
		scaleWeightsDown[6]*= qgl_weight;//nominal is with puweight
		scaleWeightsDown[7]*= qgl_weight;//nominal is with puweight
		scaleWeightsDown[8]*= qgl_weight;//nominal is with puweight
		
	
	
//		cout<<	genweight*scaleWeightsUp[current_syst]<<" "<<genweight*scaleWeightsDown[current_syst]<<endl;


			if ((current_syst==0)&&(data==1)){
				hbdtUp[current_syst]->Fill(bdt,1.);
				hbdt_atanhUp[current_syst]->Fill(TMath::ATanH((bdt+1)/2),1.);
			}else if ((current_syst==0) && (file_tag.CompareTo("interference")!=0)){
				hbdtUp[current_syst]->Fill(bdt,genweight*scaleWeightsUp[current_syst]);
				hbdt_atanhUp[current_syst]->Fill(TMath::ATanH((bdt+1)/2),genweight*scaleWeightsUp[current_syst]);
			} else if ((current_syst==0) && (file_tag.CompareTo("interference")==0)){
			///Run I		float interference_weight= 12.7733 + 1773.74/Mqq - 151127/(Mqq*Mqq) + 4.04978e+06/(Mqq*Mqq*Mqq) - 0.00044359*Mqq - 88.2666/TMath::Log(Mqq) - 1 ;
					float interference_weight= interference_kfactor * interference_func->Eval(TMath::Log(Mqq));  //Run II from Anton
					hbdtUp[current_syst]->Fill(bdt,genweight*scaleWeightsUp[current_syst]*interference_weight);
					hbdt_atanhUp[current_syst]->Fill(TMath::ATanH((bdt+1)/2),genweight*scaleWeightsUp[current_syst]*interference_weight);
			}
			if ((current_syst!=0) && (uncertainty_name[current_syst].CompareTo("JER")!=0) &&  (uncertainty_name[current_syst].CompareTo("JES")!=0) && (uncertainty_name[current_syst].CompareTo("MDG_NLO_corr")!=0   )&&(uncertainty_name[current_syst].CompareTo("int_shape")!=0)){
				hbdtUp[current_syst]->Fill(bdt,genweight*scaleWeightsUp[current_syst]);
				hbdt_atanhUp[current_syst]->Fill(TMath::ATanH((bdt+1)/2),genweight*scaleWeightsUp[current_syst]);
				hbdtDown[current_syst]->Fill(bdt,genweight*scaleWeightsDown[current_syst]);
				hbdt_atanhDown[current_syst]->Fill(TMath::ATanH((bdt+1)/2),genweight*scaleWeightsDown[current_syst]);
			} else  if ((current_syst!=0) && (uncertainty_name[current_syst].CompareTo("JES")==0)  ){
				if (passSel_JESR[0]==1) hbdtUp[current_syst]->Fill(bdt_JESR[0],genweight*scaleWeightsUp[current_syst]);
				if (passSel_JESR[0]==1) hbdt_atanhUp[current_syst]->Fill(TMath::ATanH((bdt_JESR[0]+1)/2),genweight*scaleWeightsUp[current_syst]);
				if (passSel_JESR[1]==1) hbdtDown[current_syst]->Fill(bdt_JESR[1],genweight*scaleWeightsDown[current_syst]);
				if (passSel_JESR[1]==1) hbdt_atanhDown[current_syst]->Fill(TMath::ATanH((bdt_JESR[1]+1)/2),genweight*scaleWeightsDown[current_syst]);
        }else  if ((current_syst!=0) && (uncertainty_name[current_syst].CompareTo("JER")==0)  ){
				if (passSel_JESR[2]==1) hbdtUp[current_syst]->Fill(bdt_JESR[2],genweight*scaleWeightsUp[current_syst]);
				if (passSel_JESR[2]==1) hbdt_atanhUp[current_syst]->Fill(TMath::ATanH((bdt_JESR[2]+1)/2),genweight*scaleWeightsUp[current_syst]);
				if (passSel_JESR[3]==1) hbdtDown[current_syst]->Fill(bdt_JESR[3],genweight*scaleWeightsDown[current_syst]);
				if (passSel_JESR[3]==1) hbdt_atanhDown[current_syst]->Fill(TMath::ATanH((bdt_JESR[3]+1)/2),genweight*scaleWeightsDown[current_syst]);
        } else  if ((current_syst!=0) && (uncertainty_name[current_syst].CompareTo("MDG_NLO_corr")==0)  ){
				hbdtUp[current_syst]->Fill(bdt,genweight*scaleWeightsUp[current_syst] * alpha_Mqq_up);
				hbdt_atanhUp[current_syst]->Fill(TMath::ATanH((bdt+1)/2),genweight*scaleWeightsUp[current_syst] * alpha_Mqq_up );
				hbdtDown[current_syst]->Fill(bdt,genweight*scaleWeightsDown[current_syst]*alpha_Mqq_down);
				hbdt_atanhDown[current_syst]->Fill(TMath::ATanH((bdt+1)/2),genweight*scaleWeightsDown[current_syst] * alpha_Mqq_down);
        } else  if ((current_syst!=0) && (uncertainty_name[current_syst].CompareTo("int_shape")==0)  ){
				float interference_weight_up=  interference_func_fixed->Eval(TMath::Log(Mqq));  
				float interference_weight_down= interference_down->Eval(TMath::Log(Mqq));  
			//	cout<< TMath::Log(Mqq)<< "  , "<<Mqq<<" , " << interference_weight_up<<",  "<<interference_weight_down <<endl;
		//		cout<< scaleWeightsUp[current_syst]<<" , "<<interference_weight_up<<endl;
				hbdtUp[current_syst]->Fill(bdt,genweight*scaleWeightsUp[current_syst] * interference_weight_up);
				hbdt_atanhUp[current_syst]->Fill(TMath::ATanH((bdt+1)/2),genweight*scaleWeightsUp[current_syst] * interference_weight_up );
				hbdtDown[current_syst]->Fill(bdt,genweight*scaleWeightsDown[current_syst]*interference_weight_down);
				hbdt_atanhDown[current_syst]->Fill(TMath::ATanH((bdt+1)/2),genweight*scaleWeightsDown[current_syst] *interference_weight_down );
			}

		/*	if ((current_syst==0) &&( (file_tag.CompareTo("DYJetstoLL_amc_0J")==0) || (file_tag.CompareTo("DYJetstoLL_amc_1J")==0 || (file_tag.CompareTo("DYJetstoLL_amc_2J")==0))))  {
				for (int j=0;j<30;j++){
					cout<< j<<"   "<<LHE_weights_pdf_eigen[2*j]<<"   "<<LHE_weights_pdf_eigen[2*j+1] <<endl;
					hbdt_atanhUp_pdf[j]->Fill(TMath::ATanH((bdt+1)/2),genweight*scaleWeightsUp[current_syst]*LHE_weights_pdf_eigen[2*j]);
					hbdt_atanhDown_pdf[j]->Fill(TMath::ATanH((bdt+1)/2),genweight*scaleWeightsUp[current_syst]*LHE_weights_pdf_eigen[2*j+1]);
				}
			}*/




	}

    
        hbdtUp[current_syst]->Draw();
        hbdtUp[current_syst]->Write();
        if (current_syst!=0) hbdtDown[current_syst]->Draw();
        if (current_syst!=0) hbdtDown[current_syst]->Write();
        hbdt_atanhUp[current_syst]->Draw();
        hbdt_atanhUp[current_syst]->Write();
       if (current_syst!=0)  hbdt_atanhDown[current_syst]->Draw();
       if (current_syst!=0)  hbdt_atanhDown[current_syst]->Write();
	/*	 if ((current_syst==0) &&( (file_tag.CompareTo("DYJetstoLL_amc_0J")==0) || (file_tag.CompareTo("DYJetstoLL_amc_1J")==0 || (file_tag.CompareTo("DYJetstoLL_amc_2J")==0))))  {
			for (int j=0;j<30;j++){
       		hbdt_atanhUp_pdf[j]->Draw();
      	   hbdt_atanhUp_pdf[j]->Write();
     		   hbdt_atanhDown_pdf[j]->Draw();
        	   hbdt_atanhDown_pdf[j]->Write();
			}
		}*/
}

file.Write();
file.Close();

return 0;
    
}
