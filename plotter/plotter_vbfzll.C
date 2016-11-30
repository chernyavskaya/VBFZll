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
int applyQCDScaleWeight = atoi(argv[5]);
TString QCDScaleWeight_str = std::string(argv[6]);
int applyJESWeight = atoi(argv[7]);
TString JESWeight_str = std::string(argv[8]);
TString heppyVersion = std::string(argv[9]);
TString postfix = std::string(argv[10]);
TString output = std::string(argv[11]);


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
xsec["DYJetstoLL_HT100_200"] = 173.96106; 
xsec["DYJetstoLL_HT200_400"] = 48.27802 ; 
xsec["DYJetstoLL_HT400_600"] =6.68755 ; 
xsec["DYJetstoLL_HT600_Inf"] = 2.588804; 
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

 int counter=0;


int whichQCDScaleWeight;
if ((QCDScaleWeight_str.CompareTo("none")==0)||(QCDScaleWeight_str.CompareTo("nom")==0)) whichQCDScaleWeight=0;
if (QCDScaleWeight_str.CompareTo("up")==0) whichQCDScaleWeight=1;
if (QCDScaleWeight_str.CompareTo("down")==0) whichQCDScaleWeight=2;
int whichJESWeight;
if ((JESWeight_str.CompareTo("none")==0)||(JESWeight_str.CompareTo("nom")==0)) whichJESWeight=0;
if (JESWeight_str.CompareTo("up")==0) whichJESWeight=1;
if (JESWeight_str.CompareTo("down")==0) whichJESWeight=2;

    
float gen_pos=0; 
float gen_neg=0; 
float gen_pos_weight=0; 
float gen_neg_weight=0; 


	
	Float_t presel=0;
	Float_t presel_vtype[10] = {0,0,0,0,0,0,0,0,0};
	Float_t pos_puweight=0;
	Float_t all_puweight=0.;
	Float_t puweight;
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
	Int_t events_generated;
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
 	//	events_generated = countWeighted->GetEntries();
 		if (whichQCDScaleWeight==0) events_generated = countPos->GetEntries() - countNeg->GetEntries();
 		else {
			events_generated = countLHEScale->GetBinContent( countLHEScale->FindBin( 3 + whichQCDScaleWeight) );
			if (events_generated==0) events_generated =  countPos->GetEntries() - countNeg->GetEntries();
		}
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

	Float_t bdt;
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



	if (data==1){
		genweight = 1.;
		bTagWeight = 1.;
		puweight=1.;
	}

 	
    TH1F *hVtype = new TH1F("hVtype","", 6,-1.,6.);
    hVtype->GetXaxis()->SetTitle("vtype");

	TH1F *hMqq = new TH1F("hMqq","",100.,0.,3000.);
	hMqq->GetXaxis()->SetTitle("m(qq) (GeV)");
	TH1F *hMqq_log = new TH1F("hMqq_log","",150.,0.,15.);
	hMqq_log->GetXaxis()->SetTitle("ln(m(qq)) (GeV)");
	TH1F *hqq_pt = new TH1F("hqq_pt","",80.,0.,800.);
	hqq_pt->GetXaxis()->SetTitle("p_{T}(qq) (GeV)");
	TH1F *hlepton1_pt = new TH1F("hlepton1_pt","",40.,0.,400.);
	hlepton1_pt->GetXaxis()->SetTitle("leading lepton p_{T} (GeV)");
	TH1F *hlepton2_pt = new TH1F("hlepton2_pt","",30.,0.,300.);
	hlepton2_pt->GetXaxis()->SetTitle("subleading lepton p_{T} (GeV)");
	TH1F *hlepton1_eta = new TH1F("hlepton1_eta","",40.,-4.,4.);
	hlepton1_eta->GetXaxis()->SetTitle("leading lepton #eta");
	TH1F *hlepton2_eta = new TH1F("hlepton2_eta","",40.,-4.,4.);
	hlepton2_eta->GetXaxis()->SetTitle("subleading lepton #eta");
   
	TH1F *hEtaQQ = new TH1F("hEtaQQ","",90,0.,9.);
	hEtaQQ->GetXaxis()->SetTitle("|#Delta#eta_{qq}|");
	
	TH1F *hbdt = new TH1F("hbdt","",100,-1.,1.);
	hbdt->GetXaxis()->SetTitle("BDT output");
	TH1F *hbdt_atanh = new TH1F("hbdt_atanh","",500,0.,5.);
	hbdt_atanh->GetXaxis()->SetTitle("AThanH((BDT+1)/2)");

	TH1F *hPhiQQ = new TH1F("hPhiQQ","",32,0.,3.2);
	hPhiQQ->GetXaxis()->SetTitle("|#Delta#phi_{qq}|");
    
	
	TH1F *hEtaSoftJets = new TH1F("hEtaSoftJets","",12,-3.,3.);
	hEtaSoftJets->GetXaxis()->SetTitle("|#eta^{soft}|");
    
	
	TH1F *hMassSoftJets = new TH1F("hMassSoftJets","",10,0.,100.);
	hMassSoftJets->GetXaxis()->SetTitle("m^{soft}");

	TH1F *hHTsoft = new TH1F("hHTsoft","",30,0.,300.);
	hHTsoft->GetXaxis()->SetTitle("H_{T}^{soft} (GeV)" );
	TH1F *hSoft_n2 = new TH1F("hSoft_n2","",25,0.,25.);
	hSoft_n2->GetXaxis()->SetTitle("N soft jets, p_{T} > 2 GeV");
	TH1F *hSoft_n5 = new TH1F("hSoft_n5","",10,0.,10.);
	hSoft_n5->GetXaxis()->SetTitle("N soft jets, p_{T} > 5 GeV");
	TH1F *hSoft_n10 = new TH1F("hSoft_n10","",6,0.,6.);
	hSoft_n10->GetXaxis()->SetTitle("N soft jets, p_{T} > 10 GeV");
	
	TH1F* hqgl = new TH1F("hqgl","",20.,0.,1.);
	hqgl->GetXaxis()->SetTitle("QGL 1^{st} q-jet");
	
	TH1F* hqgl2 = new TH1F("hqgl2","",20.,0.,1.);
	hqgl2->GetXaxis()->SetTitle("QGL 2^{nd} q-jet");

	
	TH1F *hPtSoftJets = new TH1F("hPtSoftJets","",30,0.,300);
	hPtSoftJets->GetXaxis()->SetTitle("p_{T}^{soft} (GeV)");
    TH1F *hPtSoftJets2 = new TH1F("hPtSoftJets2", "", 20, 0., 200.);
    hPtSoftJets2->GetXaxis()->SetTitle("2nd Soft Jet p_{T} (GeV)");
    TH1F *hPtSoftJets3 = new TH1F("hPtSoftJets3", "", 20, 0., 200.);
    hPtSoftJets3->GetXaxis()->SetTitle("3rd Soft Jet p_{T} (GeV)");
	
	TH1F *hcosOqqbb = new TH1F("hcosOqqbb","",100,-1.,1.);
	hcosOqqbb->GetXaxis()->SetTitle("cos(#theta_{bb_qq})");
	TH1F *hEtaQB1 = new TH1F("hEtaQB1","",160.,-8.,8.);
	hEtaQB1->GetXaxis()->SetTitle("#Delta#eta_{qb}^{forward}");
	TH1F *hEtaQB2 = new TH1F("hEtaQB2","",160.,-8.,8.);
	hEtaQB2->GetXaxis()->SetTitle("#Delta#eta_{qb}^{backward}");
	TH1F *hPhiQB1 = new TH1F("hPhiQB1","",32,0.,3.2);
	hPhiQB1->GetXaxis()->SetTitle("#Delta#phi_{qb}^{forward}");
	TH1F *hPhiQB2 = new TH1F("hPhiQB2","",32,0.,3.2);
	hPhiQB2->GetXaxis()->SetTitle("#Delta#phi_{qb}^{backward}");
	TH1F *hx1 = new TH1F("hx1","",100.,0.,1.);
	hx1->GetXaxis()->SetTitle("x_{1}");
	TH1F *hx2 = new TH1F("hx2","",100.,0.,1.);
	hx2->GetXaxis()->SetTitle("x_{2}");
	TH1F *hVB1_mass = new TH1F("hVB1_mass","",100,0.,1000.);
	hVB1_mass->GetXaxis()->SetTitle("M_{W'_{1}} (GeV)");
	TH1F *hVB2_mass = new TH1F("hVB2_mass","",100.,0.,1000.);
	hVB2_mass->GetXaxis()->SetTitle("M_{W'_{2}} (GeV)");

	TH1F* hEtot = new TH1F("hEtot","",150.,0.,6000.);
	hEtot->GetXaxis()->SetTitle("E^{tot} (GeV)");
	TH1F* hPxtot= new TH1F("hPxtot","",100,-500.,500.);
	hPxtot->GetXaxis()->SetTitle("P_{x}^{tot} (GeV)");
	TH1F* hPytot= new TH1F("hPytot","",100,-500.,500.);
	hPytot->GetXaxis()->SetTitle("P_{y}^{tot} (GeV)");
	TH1F* hPztot= new TH1F("hPztot","",100,-5000.,5000);
	hPztot->GetXaxis()->SetTitle("P_{z}^{tot} (GeV)");

	
	TH1F *hPtqqll = new TH1F("hPtqqll","",50.,0.,500.);
	hPtqqll->GetXaxis()->SetTitle("p_{T} of qqll system (GeV)");
	TH1F *hPhiqqll = new TH1F("hPhiqqll","",32,-3.2,3.2);
	hPhiqqll->GetXaxis()->SetTitle("-#phi of qqll system");
	TH1F *hEtaqqll = new TH1F("hEtaqqll","",160,0,8);
	hEtaqqll->GetXaxis()->SetTitle("#eta of qqll system");

	TH1F *hnPVs = new TH1F("hPVs","",50,0,50);
	hnPVs->GetXaxis()->SetTitle("nPVs");
	

	TH1F* hV_mass = new TH1F("hV_mass","",20,70.,110.);
	hV_mass->GetXaxis()->SetTitle("V_mass (GeV)");

	TH1F* hZll_mass = new TH1F("hZll_mass","",20,70.,110.);
	hZll_mass->GetXaxis()->SetTitle("m(Z) (GeV)");
	TH1F* hZll_pt = new TH1F("hZll_pt","",40,0.,400.);
	hZll_pt->GetXaxis()->SetTitle("p_{T}(Z) (GeV)");
	TH1F* hZll_eta = new TH1F("hZll_eta","",20,-5,5.);
	hZll_eta->GetXaxis()->SetTitle("#eta(Z)");
	TH1F* hZll_phi = new TH1F("hZll_phi","",32,-3.2,3.2);
	hZll_phi->GetXaxis()->SetTitle("#phi(Z)");

	TH1F *hJet1q_pt = new TH1F("hJet1q_pt","",235,30,500.);
	hJet1q_pt->GetXaxis()->SetTitle("p_{T} 1^{st} q-jet");
	TH1F *hJet1q_eta = new TH1F("hJet1q_eta","",20,-5,5);
	hJet1q_eta->GetXaxis()->SetTitle("#eta 1^{st} q-jet");
	TH1F *hJet1q_ptd = new TH1F("hJet1q_ptd","",100,0,1);
	hJet1q_ptd->GetXaxis()->SetTitle("ptd 1^{st} q-jet");
	TH1F *hJet1q_axis2= new TH1F("hJet1q_axis2","",80,0.,0.16);
	hJet1q_axis2->GetXaxis()->SetTitle("#sigma_{2} 1^{st} q-jet");
	TH1F *hJet1q_mult= new TH1F("hJet1q_mult","",30,0,30.);
	hJet1q_mult->GetXaxis()->SetTitle("N 1^{st} q-jet");
	TH1F *hJet1q_leadTrackPt= new TH1F("hJet1q_leadTrackPt","",20,0,100.);
	hJet1q_leadTrackPt->GetXaxis()->SetTitle("leading track p_{T} 1^{st} q-jet");
	TH1F *hJets12_pt = new TH1F("hJets12_pt","",335,30,600.);
	hJets12_pt->GetXaxis()->SetTitle("|p_{T}_{1q} + p_{T}_{2q}| (GeV)");
	TH1F *hJets12_pt_log = new TH1F("hJets12_pt_log","",100,0,10);
	hJets12_pt_log->GetXaxis()->SetTitle("ln|p_{T}_{1q} + p_{T}_{2q}| (GeV)");


	TH1F *hJet2q_pt = new TH1F("hJet2q_pt","",235,30,500);
	hJet2q_pt->GetXaxis()->SetTitle("p_{T} 2^{nd} q-jet");
	TH1F *hJet2q_eta = new TH1F("hJet2q_eta","",20,-5,5);
	hJet2q_eta->GetXaxis()->SetTitle("#eta 2^{nd} q-jet");
	TH1F *hJet2q_ptd = new TH1F("hJet2q_ptd","",100,0,1);
	hJet2q_ptd->GetXaxis()->SetTitle("ptd 2^{nd} q-jet");
	TH1F *hJet2q_axis2= new TH1F("hJet2q_axis2","",80,0.,0.16);
	hJet2q_axis2->GetXaxis()->SetTitle("#sigma_{2} 2^{nd} q-jet");

	TH1F *hist_bins = new TH1F("bins","",80,0.,0.16);
	
	TH1F *hJet2q_pt_log = new TH1F("hJet2q_pt_log","",80,0,8);
	hJet2q_pt_log->GetXaxis()->SetTitle("log p_{T} 2^{nd} q-jet");
	TH1F *hJet1q_pt_log = new TH1F("hJet1q_pt_log","",80,0,8);
	hJet1q_pt_log->GetXaxis()->SetTitle("log p_{T} 1^{st} q-jet");
	


	TH1F *hJet2q_mult= new TH1F("hJet2q_mult","",30,0,30.);
	hJet2q_mult->GetXaxis()->SetTitle("N 2^{nd} q-jet");
	TH1F *hJet2q_leadTrackPt= new TH1F("hJet2q_leadTrackPt","",20,0,100.);
	hJet2q_leadTrackPt->GetXaxis()->SetTitle("leading track p_{T} 2^{nd} q-jet");
	
	TH1F *hJet3_pt = new TH1F("hJet3_pt","",17,30,200);
	hJet3_pt->GetXaxis()->SetTitle("p_{T} 3^{rd} jet");


	TH1F *hmet = new TH1F("hmet","",40,0.,400.);
	hmet->GetXaxis()->SetTitle("MET p_{T} (GeV)");
	TH1F *hrho = new TH1F("hrho","",60,0.,30.);
	hrho->GetXaxis()->SetTitle("rho");
	TH1F *hHT = new TH1F("hHT","",50,0.,1000.);
	hHT->GetXaxis()->SetTitle("lhe H_{T} (GeV)" );
	TH1F *hlheHT_log = new TH1F("hlheHT_log","",100,0.,10.);
	hlheHT_log->GetXaxis()->SetTitle("ln(lhe H_{T}) (GeV)" );


	TH1F *hDeltaRelQQ = new TH1F("hDeltaRelQQ","",25.,0.,1.);
	hDeltaRelQQ->GetXaxis()->SetTitle("#Delta_{rel}(qq)");
	TH1F *hRptHard = new TH1F("hRptHard","",25.,0.,1.);
	hRptHard->GetXaxis()->SetTitle("R(p_{T}^{hard})");
	TH1F *hEtaQQSum = new TH1F("hEtaQQSum","",90,0.,9.);
	hEtaQQSum->GetXaxis()->SetTitle("|#Delta_{q_{1}}| + |#Delta_{q_{2}}| ");
	TH1F *hPhiZQ1 = new TH1F("hPhiZQ1","",32,0.,3.2);
	hPhiZQ1->GetXaxis()->SetTitle("|#Delta#phi(Z{q_{1}})|");
	TH1F* hZll_y = new TH1F("hZll_y","",20,-5,5.);
	hZll_y->GetXaxis()->SetTitle("y(Z)");
	TH1F* hZll_ystar = new TH1F("hZll_ystar","",20,-5,5.);
	hZll_ystar->GetXaxis()->SetTitle("y*(Z)");
	TH1F* hZll_zstar = new TH1F("hZll_zstar","",15,0,3.);
	hZll_zstar->GetXaxis()->SetTitle("z*(Z)");
	TH1F* hlheV_pt = new TH1F("hlheV_pt","",40,0.,400.);
	hlheV_pt->GetXaxis()->SetTitle("lheV_pt (GeV)");




   		const int numArray= 52; 
   		TH1F* histArray[numArray] = { hMqq, hEtaQQ,hHTsoft,hSoft_n2,hSoft_n5,hSoft_n10,hnPVs, hJet1q_pt, hJet1q_eta, hJet1q_ptd, hJet1q_axis2, hJet1q_mult, hJet2q_pt, hJet2q_eta, hJet2q_ptd, hJet2q_axis2, hJet2q_mult, hmet,   hJet1q_leadTrackPt, hJet2q_leadTrackPt, hqq_pt,hV_mass, hqgl, hqgl2, hZll_mass, hZll_pt, hZll_phi, hZll_eta, hrho, hlepton1_pt, hlepton2_pt, hlepton1_eta, hlepton2_eta, hHT, hDeltaRelQQ, hRptHard, hEtaQQSum, hPhiZQ1, hZll_y, hZll_ystar, hZll_zstar, hMqq_log, hlheV_pt, hJet3_pt, hlheHT_log, hPhiQQ, hJets12_pt_log, hJets12_pt, hJet1q_pt_log, hJet2q_pt_log, hbdt, hbdt_atanh };
			for (int i=0;i<numArray;i++){
				histArray[i]->Sumw2();
			}
	
		TString cut_flow_names[30] = {"triggers","2jets events","q1_pt>50","q2_pt>30","Mqq>200","leptons_pt<20","(mll-mz)<15"};
		Float_t cut_flow[30] = {0,0,0,0,0,0,0};

	float qq_matching = 0;
	float qq_matching_all = 0;
	

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


	for (int entry=0; entry<nentries;++entry){
//	for (int entry=0; entry<100;++entry){
        tree_initial->GetEntry(entry);

		for (int i=0;i<nJets;i++){
			if (whichJESWeight==1) Jet.pt[i]*= Jet.JEC_corr_up[i]/Jet.JEC_corr[i]; 	
			if (whichJESWeight==2) Jet.pt[i]*= Jet.JEC_corr_down[i]/Jet.JEC_corr[i];	
		}

		if (JSON!=1) {
			continue;
		}


		if (region.CompareTo("mu")==0) if (!(v_type==0)) continue;
		if (region.CompareTo("el")==0) if (!(v_type==1)) continue;

			

		if (data==1) PU=1.;
		else PU=puweight;
		genweight0 = genweight/TMath::Abs(genweight);
		genweight=genweight/TMath::Abs(genweight)*PU;
		genweight/=events_generated/xsec[file_tag]; 
		if  ((data!=1 ) && (file_tag.CompareTo("WW")!=0) && (file_tag.CompareTo("ZZ")!=0) && (file_tag.CompareTo("WZ")!=0)) if (whichQCDScaleWeight==1)  genweight*=LHE_weights_scale_wgt[4];
		if  ((data!=1 ) && (file_tag.CompareTo("WW")!=0) && (file_tag.CompareTo("ZZ")!=0) && (file_tag.CompareTo("WZ")!=0)) if (whichQCDScaleWeight==2)  genweight*=LHE_weights_scale_wgt[5];

		
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


		if (data==1) { 
			string file_tag_str = file_tag.Data();
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
	//	if  ( (file_tag.CompareTo("DYJetstoLL_HT100_200")==0) ||(file_tag.CompareTo("DYJetstoLL_HT200_400")==0) || (file_tag.CompareTo("DYJetstoLL_HT400_600")==0) ||(file_tag.CompareTo("DYJetstoLL_HT600_Inf")==0)   ) genweight*=ptWeightQCD(nGenVbosons, lheHT, GenVbosons_pdgId_first);
	//	if (lheHT>60) if  ((file_tag.CompareTo("DYJetstoLL")==0)  || (file_tag.CompareTo("DYJetstoLL_HT100")==0) ||(file_tag.CompareTo("DYJetstoLL_HT100_200")==0) ||(file_tag.CompareTo("DYJetstoLL_HT200_400")==0) || (file_tag.CompareTo("DYJetstoLL_HT400_600")==0) ||(file_tag.CompareTo("DYJetstoLL_HT600_Inf")==0)   ) genweight*=func_lheHT->Eval(lheHT);		
	//	if ((TMath::Log(lheHT)>4.2) && (TMath::Log(lheHT)<7.8)) if  ((file_tag.CompareTo("DYJetstoLL")==0)  || (file_tag.CompareTo("DYJetstoLL_HT100")==0) ||(file_tag.CompareTo("DYJetstoLL_HT100_200")==0) ||(file_tag.CompareTo("DYJetstoLL_HT200_400")==0) || (file_tag.CompareTo("DYJetstoLL_HT400_600")==0) ||(file_tag.CompareTo("DYJetstoLL_HT600_Inf")==0)   ) genweight*=func_lheHT->Eval(TMath::Log(lheHT));		
//		if  ((file_tag.CompareTo("DYJetstoLL")==0)  || (file_tag.CompareTo("DYJetstoLL_HT100")==0) ||(file_tag.CompareTo("DYJetstoLL_HT100_200")==0) ||(file_tag.CompareTo("DYJetstoLL_HT200_400")==0) || (file_tag.CompareTo("DYJetstoLL_HT400_600")==0) ||(file_tag.CompareTo("DYJetstoLL_HT600_Inf")==0)   ) genweight*=func_EtaQQ->Eval(qqDeltaEta);		
//
//
		Float_t Mqq_log = TMath::Log(Mqq);	
		float jets_ptSum =TMath::Log(jets_pv[0].Pt() + jets_pv[1].Pt());
		if (region.CompareTo("el")==0) {
			func_JetsPt->FixParameter(0,-467.774);
			func_JetsPt->FixParameter(1,408.471);
			func_JetsPt->FixParameter(2,-141.63);
			func_JetsPt->FixParameter(3,24.4376);
			func_JetsPt->FixParameter(4,-2.0985);
			func_JetsPt->FixParameter(5,0.0717106);

			func_EtaQQ->FixParameter(0,0.939881);
			func_EtaQQ->FixParameter(1,0.214880);
			func_EtaQQ->FixParameter(2,-0.456008);
			func_EtaQQ->FixParameter(3,0.352969);
			func_EtaQQ->FixParameter(4,-0.125972);
			func_EtaQQ->FixParameter(5,0.0229632);
			func_EtaQQ->FixParameter(6,-0.00208432);
			func_EtaQQ->FixParameter(7,0.0000746385);

	//		if ((jets_ptSum>=4.3) && (jets_ptSum<=7.7113)) if ((file_tag.CompareTo("DYJetstoLL")==0)  || (file_tag.CompareTo("DYJetstoLL_HT100")==0) ||(file_tag.CompareTo("DYJetstoLL_HT100_200")==0) ||(file_tag.CompareTo("DYJetstoLL_HT200_400")==0) || (file_tag.CompareTo("DYJetstoLL_HT400_600")==0) ||(file_tag.CompareTo("DYJetstoLL_HT600_Inf")==0)   ) genweight*=func_JetsPt->Eval(TMath::Log(jets_pv[0].Pt() + jets_pv[1].Pt()));		
	//		if (qqDeltaEta<=8.43) if ((file_tag.CompareTo("DYJetstoLL")==0)  || (file_tag.CompareTo("DYJetstoLL_HT100")==0) ||(file_tag.CompareTo("DYJetstoLL_HT100_200")==0) ||(file_tag.CompareTo("DYJetstoLL_HT200_400")==0) || (file_tag.CompareTo("DYJetstoLL_HT400_600")==0) ||(file_tag.CompareTo("DYJetstoLL_HT600_Inf")==0)   ) genweight*=func_EtaQQ->Eval(qqDeltaEta);	
			if ((Mqq_log<8.176 )&&(Mqq_log>5.256)) if ((file_tag.CompareTo("DYJetstoLL")==0)  || (file_tag.CompareTo("DYJetstoLL_HT100")==0) ||(file_tag.CompareTo("DYJetstoLL_HT100_200")==0) ||(file_tag.CompareTo("DYJetstoLL_HT200_400")==0) || (file_tag.CompareTo("DYJetstoLL_HT400_600")==0) ||(file_tag.CompareTo("DYJetstoLL_HT600_Inf")==0)   ) genweight*=func_Mqq->Eval(Mqq_log);		
		}
		if (region.CompareTo("mu")==0) {
			func_JetsPt->FixParameter(0,-328.522);
			func_JetsPt->FixParameter(1,280.199 );
			func_JetsPt->FixParameter(2,-94.6805);
			func_JetsPt->FixParameter(3,15.9021);
			func_JetsPt->FixParameter(4,-1.32777);
			func_JetsPt->FixParameter(5,0.0440635);

			func_EtaQQ->FixParameter(0,0.933843);
			func_EtaQQ->FixParameter(1,0.139400);
			func_EtaQQ->FixParameter(2,-0.336039);
			func_EtaQQ->FixParameter(3,0.297989);
			func_EtaQQ->FixParameter(4,-0.119395);
			func_EtaQQ->FixParameter(5,0.0242497);
			func_EtaQQ->FixParameter(6,-0.0024353);
			func_EtaQQ->FixParameter(7,0.000095607);

	//		if ((jets_ptSum>=4.3) && (jets_ptSum<=8.145)) if ((file_tag.CompareTo("DYJetstoLL")==0)  || (file_tag.CompareTo("DYJetstoLL_HT100")==0) ||(file_tag.CompareTo("DYJetstoLL_HT100_200")==0) ||(file_tag.CompareTo("DYJetstoLL_HT200_400")==0) || (file_tag.CompareTo("DYJetstoLL_HT400_600")==0) ||(file_tag.CompareTo("DYJetstoLL_HT600_Inf")==0)   ) genweight*=func_JetsPt->Eval(TMath::Log(jets_pv[0].Pt() + jets_pv[1].Pt()));		
		}
	//		if (qqDeltaEta<=8.21) if ((file_tag.CompareTo("DYJetstoLL")==0)  || (file_tag.CompareTo("DYJetstoLL_HT100")==0) ||(file_tag.CompareTo("DYJetstoLL_HT100_200")==0) ||(file_tag.CompareTo("DYJetstoLL_HT200_400")==0) || (file_tag.CompareTo("DYJetstoLL_HT400_600")==0) ||(file_tag.CompareTo("DYJetstoLL_HT600_Inf")==0)   ) genweight*=func_EtaQQ->Eval(qqDeltaEta);		
			if ((Mqq_log<8.176 )&&(Mqq_log>5.256)) if ((file_tag.CompareTo("DYJetstoLL")==0)  || (file_tag.CompareTo("DYJetstoLL_HT100")==0) ||(file_tag.CompareTo("DYJetstoLL_HT100_200")==0) ||(file_tag.CompareTo("DYJetstoLL_HT200_400")==0) || (file_tag.CompareTo("DYJetstoLL_HT400_600")==0) ||(file_tag.CompareTo("DYJetstoLL_HT600_Inf")==0)   ) genweight*=func_Mqq->Eval(Mqq_log);		


		float qter1 = 1.0;
		float qter2 = 1.0;
		float mu_correction1 = 1.0;
		float mu_correction2 = 1.0;


		if (data!=1) 	if (region.CompareTo("mu")==0) {
			rmcor->momcor_mc(lepton1, vLeptons_charge[0], vLeptons_trackerLayers[0], qter1);
			rmcor->momcor_mc(lepton2, vLeptons_charge[idx_2ndLepton], vLeptons_trackerLayers[idx_2ndLepton], qter2);
			}
		if (data==1) 	if (region.CompareTo("mu")==0){
			rmcor->momcor_data(lepton1, vLeptons_charge[0],  0, qter1);
			rmcor->momcor_data(lepton2, vLeptons_charge[idx_2ndLepton], 0, qter2);
			}


	
		Zll = lepton1+lepton2;	
		Float_t Zll_mass = Zll.M();
		Float_t Zll_pt = Zll.Pt();
		Float_t Zll_eta = Zll.Eta();
		Float_t Zll_phi = Zll.Phi();

		cut_flow[0]+=genweight;

		if (good_jets<2) continue;
		cut_flow[1]+=genweight;
		if (Qjet1.Pt() < 50) continue;
		cut_flow[2]+=genweight;
		if (Qjet2.Pt() < 30) continue;
		cut_flow[3]+=genweight;
		if (Mqq<200) continue;
		cut_flow[4]+=genweight;
	
	 	if (region.CompareTo("mu")==0) {
			if (lepton1.Pt()<30) continue;	
			if (lepton2.Pt()<20) continue;
		}	
	 	if (region.CompareTo("el")==0) {
			if (lepton1.Pt()<30) continue;	
			if (lepton2.Pt()<20) continue;
		}	
		cut_flow[5]+=genweight;
	 	if (region.CompareTo("el")==0) {
			if (TMath::Abs(lepton1.Eta())>2.1) continue;	
			if (TMath::Abs(lepton2.Eta())>2.1) continue;
		}	
	 	if (region.CompareTo("mu")==0) {
			if (TMath::Abs(lepton1.Eta())>2.4) continue;	
			if (TMath::Abs(lepton2.Eta())>2.4) continue;
		}	
		if (Zll_mass < 50 ) continue;
		if (TMath::Abs(Zll_mass - 91.1876)>15) continue;
		cut_flow[6]+=genweight;

		cout<<genweight<<endl;

		presel+=genweight;

		presel_vtype[(int)(v_type+1)]+=genweight;

		float Zll_ystar = Zll.Rapidity() - (Qjet1.Rapidity() + Qjet2.Rapidity()) ;
		float Zll_zstar = TMath::Abs( Zll_ystar/ (Qjet1.Rapidity()-Qjet2.Rapidity() )) ;
		Float_t DeltaEtaQQSum = TMath::Abs(Qjet1.Eta()) +  TMath::Abs(Qjet2.Eta());
		Float_t PhiZQ1 = TMath::Abs(Zll.DeltaPhi(Qjet1));
		Float_t DeltaRelQQ = (Qjet1+Qjet2).Pt()/( Qjet1.Pt()+Qjet2.Pt()); 
		Float_t RptHard = (Qjet1+Qjet2+ Zll).Pt()/( Qjet1.Pt()+Qjet2.Pt() + Zll.Pt()); 



/////////////////////////////////////////////////////////


			counter++;

            hVtype->Fill(v_type,genweight);
		   	hMqq->Fill(Mqq,genweight);
		   	hMqq_log->Fill(TMath::Log(Mqq),genweight);
			   hqq_pt->Fill(qq_pt,genweight);
			   hEtaQQ->Fill(qqDeltaEta,genweight);
		 	   hPhiQQ->Fill(qqDeltaPhi,genweight);
		   	hZll_mass->Fill(Zll_mass,genweight);
		   	hZll_pt->Fill(Zll_pt,genweight);
		   	hZll_phi->Fill(Zll_phi,genweight);
		   	hZll_eta->Fill(Zll_eta,genweight);
			   hHTsoft->Fill(Jet.HTsoft,genweight);
			   hSoft_n2->Fill(Jet.nsoft2, genweight);
			   hSoft_n5->Fill(Jet.nsoft5, genweight);
			   hSoft_n10->Fill(Jet.nsoft10, genweight);
				hnPVs->Fill(nPVs,genweight);
				hqgl->Fill(Jet.qgl[jets_indices[0]],genweight);
				hqgl2->Fill(Jet.qgl[jets_indices[1]],genweight);
			
				hJet1q_pt->Fill(jets_pv[0].Pt(),genweight);
				hJet1q_eta->Fill(jets_pv[0].Eta(),genweight);
				hJet1q_ptd->Fill(Jet.ptd[jets_indices[0]],genweight);
				hJet1q_axis2->Fill(TMath::Exp((-1)*Jet.axis2[jets_indices[0]]),genweight);
				hJet1q_mult->Fill(Jet.mult[jets_indices[0]],genweight);
				hJet1q_leadTrackPt->Fill(Jet.leadTrackPt[jets_indices[0]],genweight);
				hJet2q_pt->Fill(jets_pv[1].Pt(),genweight);
				hJet2q_eta->Fill(jets_pv[1].Eta(),genweight);
				hJet2q_ptd->Fill(Jet.ptd[jets_indices[1]],genweight);
				hJet2q_axis2->Fill(TMath::Exp((-1)*Jet.axis2[jets_indices[1]]),genweight);
				hJet2q_mult->Fill(Jet.mult[jets_indices[1]],genweight);
				hJet2q_leadTrackPt->Fill(Jet.leadTrackPt[jets_indices[1]],genweight);
				hmet->Fill(met_pt,genweight);
				hV_mass->Fill(V_mass,genweight);
				hrho->Fill(rho,genweight);
				hlepton1_pt->Fill(lepton1.Pt(),genweight);
				hlepton2_pt->Fill(lepton2.Pt(),genweight);
				hlepton1_eta->Fill(lepton1.Eta(),genweight);
				hlepton2_eta->Fill(lepton2.Eta(),genweight);
				hHT->Fill(lheHT ,genweight);
				hlheHT_log->Fill(TMath::Log(lheHT) ,genweight);
				
				hDeltaRelQQ->Fill(DeltaRelQQ,genweight);
				hRptHard->Fill(RptHard,genweight);
				hEtaQQSum->Fill(DeltaEtaQQSum,genweight);
				hPhiZQ1->Fill(PhiZQ1,genweight);
				hZll_y->Fill(Zll.Rapidity(),genweight);
				hZll_ystar->Fill(Zll_ystar   ,genweight);
				hZll_zstar->Fill(Zll_zstar,genweight);
				hlheV_pt->Fill(lheV_pt  ,genweight);
				hJet3_pt->Fill(jet3_pt ,genweight);	
			
				hJets12_pt->Fill((jets_pv[0].Pt() + jets_pv[1].Pt()),genweight);
				hJets12_pt_log->Fill(TMath::Log(jets_pv[0].Pt() + jets_pv[1].Pt()),genweight);
				hJet1q_pt_log->Fill(TMath::Log(jets_pv[0].Pt()),genweight);
				hJet2q_pt_log->Fill(TMath::Log(jets_pv[1].Pt()),genweight);

				hbdt->Fill(bdt,genweight);
				hbdt_atanh->Fill(TMath::ATanH((bdt+1)/2),genweight);

		
		if (genweight>0) pos_weight_presel++;
		float mcweight=genweight*events_generated/xsec[file_tag]; 

		if (genweight>0) gen_pos_weight+=mcweight;
		if (genweight<0) gen_neg_weight+=mcweight;
		if (genweight>0) gen_pos+=genweight0;
		if (genweight<0) gen_neg+=genweight0;

				
			global_counter++;
        }

		cout<<counter<<endl;
		TFile file(output+"/"+file_tag+"_"+region+"_QCDScale"+QCDScaleWeight_str+"_JES"+JESWeight_str+"_"+heppyVersion+"_"+postfix+".root","recreate");
    
		for (int i=0;i<numArray;++i){
    	    	histArray[i]->SetLineWidth(2);
    	   	histArray[i]->GetYaxis()->SetTitle("N_{events}");
       		histArray[i]->GetYaxis()->SetTitleFont(42);
       		histArray[i]->GetYaxis()->SetTitleSize(0.060);
        		histArray[i]->GetYaxis()->SetTitleOffset(0.8);
        		histArray[i]->SetLineColor(kBlue);
        		histArray[i]->Draw();
        		histArray[i]->Write();
   		} 
    		file.Write();
    		file.Close();
	 ofstream out(output+"/"+file_tag+"_"+region+"_QCDScale"+QCDScaleWeight_str+"_JES"+JESWeight_str+"_"+heppyVersion+"_"+postfix+".txt");
	out<< "positive pure selected = "<<gen_pos<<"  , positive weighted selected =  "<<gen_pos_weight<<" , negative pure selected = "<<gen_neg<< ", negative weighted selected = "<<gen_neg_weight<< ", all evetns in the begining = "<<events_generated<<" , xsec = "<<xsec[file_tag]<<endl;
	out<<"positive weight in so many events : "<<  pos_weight_presel<<endl;
	for (int i=0;i<7;i++)
		out<<cut_flow_names[i]<<"\t";
	out<<endl;
	for (int i=0;i<7;i++){
		cut_flow[i]=cut_flow[i];
		out<<cut_flow[i]<<"\t";
	}
	out<<endl;

return 0;
    
}
