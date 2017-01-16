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
//#include "/afs/cern.ch/work/n/nchernya/VBFZll/plotter/muon_corrections/RoccoR.cc"
//#include "/afs/cern.ch/work/n/nchernya/VBFZll/plotter/muon_corrections/RoccoR.h"
//#include "/afs/cern.ch/work/n/nchernya/VBFZll/plotter/muon_corrections/rochcor2016.cc"
//#include "/afs/cern.ch/work/n/nchernya/VBFZll/plotter/muon_corrections/rochcor2016.h"
//#include "EWcorr.C"
#include "/mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/rochcor2016.h"
#include "/mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/rochcor2016.cc"
#include "/mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/RoccoR.cc"
#include "/mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/RoccoR.h"

Double_t erf( Double_t *x, Double_t *par){
  return par[0]/2.*(1.+TMath::Erf((x[0]-par[1])/par[2]));
}
Double_t erf2( Double_t *x, Double_t *par){
  return par[0]/2.*(1.+TMath::Erf((x[0]-par[1])/par[2]))+ (1.-par[0]);
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
	Int_t id[njets];
	Int_t puId[njets];
	Float_t HTsoft;
	Float_t pt_regVBF[njets];	
	Float_t ptd[njets];
	Float_t axis2[njets];
	Int_t mult[njets];
	Float_t leadTrackPt[njets];
	Float_t blike_VBF[njets];
} Jets;

using namespace std;


int main(int argc, char* argv[]){

//gROOT->ProcessLine(".L /afs/cern.ch/work/n/nchernya/VBFZll/plotter/muon_corrections/./RoccoR_cc.so");
//gROOT->ProcessLine(".L /afs/cern.ch/work/n/nchernya/VBFZll/plotter/muon_corrections/./rochcor2016_cc.so");
//gROOT->ProcessLine(".L /afs/cern.ch/work/n/nchernya/VBFZll/plotter/muon_corrections/RoccoR.cc++");
//gROOT->ProcessLine(".L /afs/cern.ch/work/n/nchernya/VBFZll/plotter/muon_corrections/rochcor2016.cc++");
gROOT->ProcessLine(".L /mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/RoccoR.cc++");
gROOT->ProcessLine(".L /mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/rochcor2016.cc++");



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
//xsec["DYJetstoLL_HT100_200"] = 147.40 ; 
//xsec["DYJetstoLL_HT200_400"] = 40.99 ; 
//xsec["DYJetstoLL_HT400_600"] = 5.678 ; 
//xsec["DYJetstoLL_HT600_Inf"] = 2.198; 
xsec["DYJetstoLL_HT100_200"] = 181.302; 
xsec["DYJetstoLL_HT200_400"] =50.4177  ; 
xsec["DYJetstoLL_HT400_600"] =6.98394; 
xsec["DYJetstoLL_HT600_Inf"] =2.70354 ; 
xsec["DYJetstoLL_Pt-100_amc"] = 5765.4; 
xsec["DYJetstoLL_Pt-100To250_amc"] = 83.12; 
xsec["DYJetstoLL_Pt-250To400_amc"] =3.047 ; 
xsec["DYJetstoLL_Pt-400To650_amc"] = 0.3921 ; 
xsec["DYJetstoLL_Pt-650ToInf_amc"] = 0.03636 ; 
xsec["TT"] =809.;
xsec["WW"] =118.7;
xsec["WZ"] = 47.13;
xsec["ZZ"] =16.523;
xsec["EWK_LLJJ"]=1.664;
xsec["interference"]=1.664;

xsec["WJetsToLNu_amc"]  = 61526.7; //not going to use these ones
xsec["WJetsToLNu"]  = 61526.7;
xsec["WJetsToLNu_HT100"]  = 61526.7;
xsec["WJetsToLNu_HT100To200"] = 1627.45 ;
xsec["WJetsToLNu_HT200To400"] = 435.236  ;
xsec["WJetsToLNu_HT400To600"] = 59.18109;
xsec["WJetsToLNu_HT600To800"] =14.5805;
xsec["WJetsToLNu_HT800To1200"] = 6.656210;
xsec["WJetsToLNu_HT1200To2500"] = 1.608089;
xsec["WJetsToLNu_HT2500ToInf"] = 0.0389135;

xsec["ST_tW"] = 71.7 ;			//inclusive decays
xsec["ST_s-channel"] = 3.36; //leptonic decays
xsec["ST_t-channel_top_4f_inclusiveDecays"] = 136.02;
xsec["ST_t-channel_antitop_4f_inclusiveDecays"] = 80.95;
xsec["ST_t-channel_top_4f_leptonDecays"] = 44.33;  //leptonDecays  , multiplied with BR 0.325
xsec["ST_t-channel_antitop_4f_leptonDecays"] = 26.38;//leptonDecays ,multiplied with BR 0.325


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
	Int_t HLT_IsoMu27;
	Int_t HLT_IsoTkMu27;
	Int_t HLT_IsoMu24;
	Int_t HLT_Ele27_eta2p1;
	TFile *file_initial;
	TChain *tree_initial;

	Int_t nvLeptons;
	const int brLeptons=13;
	Float_t vLeptons_pt[30], vLeptons_eta[30], vLeptons_phi[30], vLeptons_mass[30], vLeptons_SF_IdCutLoose[30], vLeptons_SF_IdCutTight[30], vLeptons_SF_IsoLoose[30], vLeptons_SF_IsoTight[30],vLeptons_SF_trk_eta[30], vLeptons_SF_HLT_RunD4p2[30],vLeptons_SF_HLT_RunD4p3[30] ;
	Int_t vLeptons_charge[30], vLeptons_pdgId[30],vLeptons_trackerLayers[30]; 
	TString str_leptons[brLeptons] = {"vLeptons_pt", "vLeptons_eta", "vLeptons_phi", "vLeptons_mass", "vLeptons_charge", "vLeptons_pdgId", "vLeptons_SF_IdCutLoose", "vLeptons_SF_IdCutTight", "vLeptons_SF_IsoLoose","vLeptons_SF_IsoTight","vLeptons_SF_trk_eta","vLeptons_SF_HLT_RunD4p2","vLeptons_SF_HLT_RunD4p3"};


	
////////////////////////////
	TFile* file_trig_el = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/skimmed/TriggerEffMap_electronTriggerEfficiencyHLT_Ele27_WPLoose_eta2p1_WP90_BCDEF.root");
	TH2F* trig_el = (TH2F*)file_trig_el->Get("TriggerEffMap_electronTriggerEfficiencyHLT_Ele27_WPLoose_eta2p1_WP90_BCDEF");
	TFile* file_trig_mu_bf = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/skimmed/TriggerEffMap_MuonTrigger_data_all_IsoMu22_OR_IsoTkMu22_pteta_Run2016B_beforeL2Fix.root");
	TH2F* trig_mu_bf = (TH2F*)file_trig_mu_bf->Get("TriggerEffMap_MuonTrigger_data_all_IsoMu22_OR_IsoTkMu22_pteta_Run2016B_beforeL2Fix");
	TFile* file_trig_mu_aft = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/skimmed/TriggerEffMap_MuonTrigger_data_all_IsoMu22_OR_IsoTkMu22_pteta_Run2016B_afterL2Fix.root");
	TH2F* trig_mu_aft = (TH2F*)file_trig_mu_aft->Get("TriggerEffMap_MuonTrigger_data_all_IsoMu22_OR_IsoTkMu22_pteta_Run2016B_afterL2Fix");

	
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
 		events_generated = countWeighted->GetEntries();
		if (events_generated==0) events_generated =  countPos->GetEntries() - countNeg->GetEntries();
 	//	events_generated = countPos->GetEntries() - countNeg->GetEntries();
		events_generated_muF_QCDUp = countLHEScale->GetBinContent( countLHEScale->FindBin(0) );
		events_generated_muF_QCDDown = countLHEScale->GetBinContent( countLHEScale->FindBin(1) );
		events_generated_muR_QCDUp = countLHEScale->GetBinContent( countLHEScale->FindBin(2) );
		events_generated_muR_QCDDown = countLHEScale->GetBinContent( countLHEScale->FindBin(3) );
		events_generated_muFR_QCDUp = countLHEScale->GetBinContent( countLHEScale->FindBin(4) );
		events_generated_muFR_QCDDown = countLHEScale->GetBinContent( countLHEScale->FindBin(5) );
	} else events_generated = 1;
    Jets Jet;
    Float_t v_type;
    Float_t wrong_type=0.;
    Int_t nJets;
	Float_t JSON;	
	Float_t nPVs;	
	Float_t rho;	
	Float_t lheHT;	
	Float_t lheV_pt;	

	const int Nsyst = 8;
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
Float_t LHE_weights_scale_wgt[10];
	
	float V_mass;

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
	tree_initial->SetBranchAddress("met_pt",&met_pt);
	tree_initial->SetBranchAddress("met_phi",&met_phi);
	
	tree_initial->SetBranchAddress("softActivityJets_pt",Jet.soft_pt);
	tree_initial->SetBranchAddress("softActivityJets_eta",Jet.soft_eta);
	tree_initial->SetBranchAddress("softActivityJets_mass",Jet.soft_mass);
	tree_initial->SetBranchAddress("softActivity_HT",&Jet.HTsoft);
	tree_initial->SetBranchAddress("softActivity_njets2",&Jet.nsoft2);
	tree_initial->SetBranchAddress("softActivity_njets5",&Jet.nsoft5);
	tree_initial->SetBranchAddress("softActivity_njets10",&Jet.nsoft10);
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
 	tree_initial->SetBranchAddress("HLT_BIT_HLT_IsoMu27_v",&HLT_IsoMu27);
 	tree_initial->SetBranchAddress("HLT_BIT_HLT_IsoMu24_v",&HLT_IsoMu24);
 	tree_initial->SetBranchAddress("HLT_BIT_HLT_IsoTkMu27_v",&HLT_IsoTkMu27);
 	tree_initial->SetBranchAddress("HLT_BIT_HLT_Ele27_eta2p1_WPLoose_Gsf_v",&HLT_Ele27_eta2p1);
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


//	for (int i=0;i<brLeptons;i++){
//		tree_initial->SetBranchAddress(str_leptons[i],);
//	}	

	tree_initial->SetBranchAddress("LHE_weights_pdf_wgt",LHE_weights_pdf_wgt);
	tree_initial->SetBranchAddress("LHE_weights_scale_wgt",LHE_weights_scale_wgt);
	tree_initial->SetBranchAddress("lheHT",&lheHT);
	tree_initial->SetBranchAddress("lheV_pt",&lheV_pt);
	tree_initial->SetBranchAddress("BDT_VBF",&bdt);
	tree_initial->SetBranchAddress("BDT_VBF_JES_up",bdt_JESR[0]);
	tree_initial->SetBranchAddress("BDT_VBF_JES_down",bdt_JESR[1]);
	tree_initial->SetBranchAddress("BDT_VBF_JER_up",bdt_JESR[2]);
	tree_initial->SetBranchAddress("BDT_VBF_JER_down",bdt_JESR[3]);

	tree_initial->SetBranchAddress("PassSelection_JES_up",passSel_JESR[0]);
	tree_initial->SetBranchAddress("PassSelection_JES_down",passSel_JESR[1]);
	tree_initial->SetBranchAddress("PassSelection_JER_up",passSel_JESR[2]);
	tree_initial->SetBranchAddress("PassSelection_JER_down",passSel_JESR[3]);
	tree_initial->SetBranchAddress("PassSelection_nom",&passSel);



	if (data==1){
		genweight = 1.;
		bTagWeight = 1.;
		puweight=1.;
	}

	float lumi = 22000;


	TString uncertainty_name[Nsyst] = {"","puWeight","LHE_weights_scale","LHE_weights_scale_muF","LHE_weights_scale_muR","JES","JER","MDG_NLO_corr"};
	TH1F *hbdtUp[20];
	TH1F *hbdt_atanhUp[20];
	TH1F *hbdtDown[20];
	TH1F *hbdt_atanhDown[20];
	for (int i=0;i<Nsyst;i++){
		char histTitleUp[50];
		char histTitleDown[50];
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



	

	int nentries = tree_initial->GetEntries() ;
	

//	TF1 *func_lheHT = new TF1("func_lheHT","([0]+[1]*x+[2]*x*x+[3]*x*x*x)*TMath::Exp(-1*[4]*x)",60,4000);
//	func_lheHT->FixParameter(0,  -1.71063e+00);
//	func_lheHT->FixParameter(1, 6.90159e-02 );
///	func_lheHT->FixParameter(2, -2.83168e-04);
//	func_lheHT->FixParameter(3, 4.69007e-07);
//	func_lheHT->FixParameter(4, 7.07950e-03 );

	TF1* func_lheHT = new TF1("func_lheHT","pol6",4.2,7.8);
	func_lheHT->FixParameter(0,89.0139);
	func_lheHT->FixParameter(1,-275.535);
	func_lheHT->FixParameter(2,195.308);
	func_lheHT->FixParameter(3,-61.2467);
	func_lheHT->FixParameter(4,9.8217);
	func_lheHT->FixParameter(5,-0.791744);
	func_lheHT->FixParameter(6,0.0255211);
//	TF1* func_EtaQQ = new TF1("func_EtaQQ","[0]+[1]*x+[2]*x*x+[3]*x*x*x",0.,16.);
//	func_EtaQQ->FixParameter(0,1.33558e+00);
//	func_EtaQQ->FixParameter(1,-2.60070e-01 );
//	func_EtaQQ->FixParameter(2,5.87759e-02);
//	func_EtaQQ->FixParameter(3,-4.64109e-03 );

	TF1* func_JetsPt = new TF1("func_JetsPt","pol5",4.3,10);
	TF1* func_EtaQQ = new TF1("func_EtaQQ","pol7",0,10);
	TF1* func_Mqq = new TF1("func_Mqq","pol6",0,10);
	func_Mqq->FixParameter(0,5968.223851);
	func_Mqq->FixParameter(1,-5554.340558);
	func_Mqq->FixParameter(2,2146.273308);
	func_Mqq->FixParameter(3,-440.733829);
	func_Mqq->FixParameter(4,50.729466);
	func_Mqq->FixParameter(5,-3.103351);
	func_Mqq->FixParameter(6,0.0788278);

	rochcor2016 *rmcor = new rochcor2016();



TFile file(output+"/"+file_tag+"_"+region+"_"+heppyVersion+"_"+postfix+".root","recreate");
float scaleWeightsUp[20];
float scaleWeightsDown[20];
string file_tag_str = file_tag.Data();
TString uncNotAppl[20] = {"","data","WW_WZ_ZZ","WW_WZ_ZZ","WW_WZ_ZZ","","",""};
for (int current_syst=0;current_syst<Nsyst;current_syst++){
	if ((data==1)&&(current_syst!=0)) continue;
	string uncNotAppl_str = uncNotAppl[current_syst].Data();
	if (current_syst!=0) if (uncNotAppl_str.find(file_tag_str)!=std::string::npos) continue;
	if  ( (uncertainty_name[current_syst].CompareTo("MDG_NLO_corr")==0) &&(  !((file_tag.CompareTo("DYJetstoLL")==0)  || (file_tag.CompareTo("DYJetstoLL_HT100")==0) ||(file_tag.CompareTo("DYJetstoLL_HT100_200")==0) ||(file_tag.CompareTo("DYJetstoLL_HT200_400")==0) || (file_tag.CompareTo("DYJetstoLL_HT400_600")==0) ||(file_tag.CompareTo("DYJetstoLL_HT600_Inf")==0)) ) )  continue;
	bdt = 0;
	passSel=0;
	genweight = 0 ;
	for (int entry=0; entry<nentries;++entry){
        tree_initial->GetEntry(entry);

		
	
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

		if (JSON!=1) {
			continue;
		}


		if (region.CompareTo("mu")==0) if (!(v_type==0)) continue;
		if (region.CompareTo("el")==0) if (!(v_type==1)) continue;

			

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
		lepton1.SetPtEtaPhiM(vLeptons_pt[0], vLeptons_eta[0], vLeptons_phi[0], vLeptons_mass[0]);	
		int idx_2ndLepton = 0;
		for (int i=1; i<nvLeptons;i++ ){
			if (vLeptons_charge[0]*vLeptons_charge[i] < 0) {
				idx_2ndLepton=i;
				break;
			}
		}
		lepton2.SetPtEtaPhiM(vLeptons_pt[idx_2ndLepton], vLeptons_eta[idx_2ndLepton], vLeptons_phi[idx_2ndLepton], vLeptons_mass[idx_2ndLepton]);
		if (data!=1) 	if (region.CompareTo("mu")==0) {
			rmcor->momcor_mc(lepton1, vLeptons_charge[0], vLeptons_trackerLayers[0], qter1);
			rmcor->momcor_mc(lepton2, vLeptons_charge[idx_2ndLepton], vLeptons_trackerLayers[idx_2ndLepton], qter2);
			}
		if (data==1) 	if (region.CompareTo("mu")==0){
			rmcor->momcor_data(lepton1, vLeptons_charge[0],  0, qter1);
			rmcor->momcor_data(lepton2, vLeptons_charge[idx_2ndLepton], 0, qter2);
			}


		if (data==1) { 
			if  (file_tag_str.find("SingleMuon")!=std::string::npos) if (!((HLT_IsoMu27==1) || (HLT_IsoTkMu27==1)  )) continue; 
			if  (file_tag_str.find("SingleElectron")!=std::string::npos) if (!(HLT_Ele27_eta2p1 == 1)) continue;
		} else if (data!=1) {
			if (region.CompareTo("mu")==0) {
			//	cout<<vLeptons_SF_HLT_RunD4p2[0]<<"  "<<vLeptons_SF_HLT_RunD4p3[0]<<endl;
			//	genweight*=(0.032*vLeptons_SF_HLT_RunD4p2[0] + 0.96799*vLeptons_SF_HLT_RunD4p3[0]);
			//	cout<<genweight<<endl;
				float SF_mu_bf_err1 = 0.;
				float SF_mu_bf_err2 = 0.;
				float SF_mu_aft_err1 = 0.;
				float SF_mu_aft_err2 = 0.;
				bool abs=1;
				float eff1 =0.02772*getScaleFactor(trig_mu_bf, lepton1.Pt(), lepton1.Eta(), SF_mu_bf_err1,abs ) + 0.97227*getScaleFactor(trig_mu_aft, lepton1.Pt(), lepton1.Eta(), SF_mu_bf_err1,abs ) ;  	
				float eff2 =0.02772*getScaleFactor(trig_mu_bf, lepton2.Pt(), lepton2.Eta(), SF_mu_bf_err2,abs ) + 0.97227*getScaleFactor(trig_mu_aft, lepton2.Pt(), lepton2.Eta(), SF_mu_bf_err2,abs  ) ;  
				genweight*= eff1*(1-eff2)*eff1 + eff2*(1-eff1)*eff2 + eff1*eff1*eff2*eff2; 	

	
				genweight*= vLeptons_SF_IdCutLoose[0]*vLeptons_SF_IdCutLoose[1] * vLeptons_SF_IsoLoose[0]* vLeptons_SF_IsoLoose[1]* vLeptons_SF_trk_eta[0]*vLeptons_SF_trk_eta[1];
			}
			if (region.CompareTo("el")==0) {
				float SF_el_err1 = 0.;
				float SF_el_err2 = 0.;
				bool abs=0;
				float eff1 = getScaleFactor(trig_el, lepton1.Pt(), lepton1.Eta(), SF_el_err1,abs );  	
				float eff2 = getScaleFactor(trig_el, lepton2.Pt(), lepton2.Eta(), SF_el_err2,abs ); 
				genweight*= eff1*(1-eff2)*eff1 + eff2*(1-eff1)*eff2 + eff1*eff1*eff2*eff2; 	
			}
		}

		if  ((file_tag.CompareTo("DYJetstoLL")==0) 
|| (file_tag.CompareTo("DYJetstoLL_amc")==0) || (file_tag.CompareTo("DYJetstoLL_Pt-100_amc")==0)|| (file_tag.CompareTo("DYJetstoLL_Pt-100To250_amc")==0)|| (file_tag.CompareTo("DYJetstoLL_Pt-250To400_amc")==0) || (file_tag.CompareTo("DYJetstoLL_Pt-400To650_amc")==0)|| (file_tag.CompareTo("DYJetstoLL_Pt-650ToInf_amc")==0)  
 || (file_tag.CompareTo("DYJetstoLL_HT100")==0) ||(file_tag.CompareTo("DYJetstoLL_HT100_200")==0) ||(file_tag.CompareTo("DYJetstoLL_HT200_400")==0) || (file_tag.CompareTo("DYJetstoLL_HT400_600")==0) ||(file_tag.CompareTo("DYJetstoLL_HT600_Inf")==0)   ) genweight*=ptWeightEWK(nGenVbosons, GenVbosons_pt_first, VtypeSim, GenVbosons_pdgId_first); 
//
//
		Float_t Mqq_log = TMath::Log(Mqq);	
		float jets_ptSum =TMath::Log(jets_pv[0].Pt() + jets_pv[1].Pt());
		float alpha=0;
		if (region.CompareTo("el")==0) {

			if ((Mqq_log<8.176 )&&(Mqq_log>5.256)) if ((file_tag.CompareTo("DYJetstoLL")==0)  || (file_tag.CompareTo("DYJetstoLL_HT100")==0) ||(file_tag.CompareTo("DYJetstoLL_HT100_200")==0) ||(file_tag.CompareTo("DYJetstoLL_HT200_400")==0) || (file_tag.CompareTo("DYJetstoLL_HT400_600")==0) ||(file_tag.CompareTo("DYJetstoLL_HT600_Inf")==0)   ) {
				genweight*=func_Mqq->Eval(Mqq_log);	
				alpha=	func_Mqq->Eval(Mqq_log) - 1;
			}
		}
		if (region.CompareTo("mu")==0) {

			if ((Mqq_log<8.176 )&&(Mqq_log>5.256)) if ((file_tag.CompareTo("DYJetstoLL")==0)  || (file_tag.CompareTo("DYJetstoLL_HT100")==0) ||(file_tag.CompareTo("DYJetstoLL_HT100_200")==0) ||(file_tag.CompareTo("DYJetstoLL_HT200_400")==0) || (file_tag.CompareTo("DYJetstoLL_HT400_600")==0) ||(file_tag.CompareTo("DYJetstoLL_HT600_Inf")==0)   ){
				 genweight*=func_Mqq->Eval(Mqq_log);		
			    alpha=	func_Mqq->Eval(Mqq_log) - 1;
			}
		}


	
	
//		cout<<	genweight*scaleWeightsUp[current_syst]<<" "<<genweight*scaleWeightsDown[current_syst]<<endl;

			if ((current_syst==0)&&(data==1)){
				hbdtUp[current_syst]->Fill(bdt,1.);
				hbdt_atanhUp[current_syst]->Fill(TMath::ATanH((bdt+1)/2),1.);
			}else if ((current_syst==0) && (file_tag.CompareTo("interference")!=0)){
				hbdtUp[current_syst]->Fill(bdt,genweight*scaleWeightsUp[current_syst]);
				hbdt_atanhUp[current_syst]->Fill(TMath::ATanH((bdt+1)/2),genweight*scaleWeightsUp[current_syst]);
			} else if ((current_syst==0) && (file_tag.CompareTo("interference")==0)){
					float interference_weight= 12.7733 + 1773.74/Mqq - 151127/(Mqq*Mqq) + 4.04978e+06/(Mqq*Mqq*Mqq) - 0.00044359*Mqq - 88.2666/TMath::Log(Mqq) - 1 ;
					hbdtUp[current_syst]->Fill(bdt,genweight*scaleWeightsUp[current_syst]*interference_weight);
					hbdt_atanhUp[current_syst]->Fill(TMath::ATanH((bdt+1)/2),genweight*scaleWeightsUp[current_syst]*interference_weight);
			}
			if ((current_syst!=0) && (uncertainty_name[current_syst].CompareTo("JER")!=0) &&  (uncertainty_name[current_syst].CompareTo("JES")!=0) && (uncertainty_name[current_syst].CompareTo("MDG_NLO_corr")!=0   )){
				hbdtUp[current_syst]->Fill(bdt,genweight*scaleWeightsUp[current_syst]);
				hbdt_atanhUp[current_syst]->Fill(TMath::ATanH((bdt+1)/2),genweight*scaleWeightsUp[current_syst]);
				hbdtDown[current_syst]->Fill(bdt,genweight*scaleWeightsDown[current_syst]);
				hbdt_atanhDown[current_syst]->Fill(TMath::ATanH((bdt+1)/2),genweight*scaleWeightsDown[current_syst]);
			} else  if ((current_syst!=0) && (uncertainty_name[current_syst].CompareTo("JES")==0)  ){
				if (passSel[0]==1) hbdtUp[current_syst]->Fill(bdt_JESR[0],genweight*scaleWeightsUp[current_syst]);
				if (passSel[0]==1) hbdt_atanhUp[current_syst]->Fill(TMath::ATanH((bdt_JESR[0]+1)/2),genweight*scaleWeightsUp[current_syst]);
				if (passSel[1]==1) hbdtDown[current_syst]->Fill(bdt_JESR[1],genweight*scaleWeightsDown[current_syst]);
				if (passSel[1]==1) hbdt_atanhDown[current_syst]->Fill(TMath::ATanH((bdt_JESR[1]+1)/2),genweight*scaleWeightsDown[current_syst]);
        }else  if ((current_syst!=0) && (uncertainty_name[current_syst].CompareTo("JER")==0)  ){
				if (passSel[2]==1) hbdtUp[current_syst]->Fill(bdt_JESR[2],genweight*scaleWeightsUp[current_syst]);
				if (passSel[2]==1) hbdt_atanhUp[current_syst]->Fill(TMath::ATanH((bdt_JESR[2]+1)/2),genweight*scaleWeightsUp[current_syst]);
				if (passSel[3]==1) hbdtDown[current_syst]->Fill(bdt_JESR[3],genweight*scaleWeightsDown[current_syst]);
				if (passSel[3]==1) hbdt_atanhDown[current_syst]->Fill(TMath::ATanH((bdt_JESR[3]+1)/2),genweight*scaleWeightsDown[current_syst]);
        } else  if ((current_syst!=0) && (uncertainty_name[current_syst].CompareTo("MDG_NLO_corr")==0)  ){
				hbdtUp[current_syst]->Fill(bdt,genweight*scaleWeightsUp[current_syst] * (1+3./2.*alpha));
				hbdt_atanhUp[current_syst]->Fill(TMath::ATanH((bdt+1)/2),genweight*scaleWeightsUp[current_syst] * (1+3./2.*alpha));
				hbdtDown[current_syst]->Fill(bdt,genweight*scaleWeightsDown[current_syst]* (1+1./2.*alpha));
				hbdt_atanhDown[current_syst]->Fill(TMath::ATanH((bdt+1)/2),genweight*scaleWeightsDown[current_syst] * (1+1./2.*alpha));
			}




    
        hbdtUp[current_syst]->Draw();
        hbdtUp[current_syst]->Write();
        if (current_syst!=0) hbdtDown[current_syst]->Draw();
        if (current_syst!=0) hbdtDown[current_syst]->Write();
        hbdt_atanhUp[current_syst]->Draw();
        hbdt_atanhUp[current_syst]->Write();
       if (current_syst!=0)  hbdt_atanhDown[current_syst]->Draw();
       if (current_syst!=0)  hbdt_atanhDown[current_syst]->Write();
}

file.Write();
file.Close();

return 0;
    
}
