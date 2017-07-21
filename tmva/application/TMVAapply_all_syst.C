#define TMVAapply_all_syst_cxx
#include "TMVAapply_all_syst.h"
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TLorentzVector.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "/afs/cern.ch/work/n/nchernya/VBFZll/plotter/EWcorr.C"
#include "/mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/mucorr/2016/RoccoR.cc"


const int njets = 30;


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




typedef struct {
	Float_t Mqq;
	Float_t qq_pt;
	Float_t DeltaEtaQQ;
	Float_t Jet3_pt;
	Float_t axis2_jet1;
	Float_t axis2_jet2;
	Float_t Jet2q_pt; 
	Float_t Jet2q_leadTrackPt;
	Float_t Jet1q_pt;
	Float_t Jet1q_leadTrackPt;
	Float_t Zll_zstar;
	Float_t RptHard;
	Float_t Zll_pt;
}TMVAstruct;

float getScaleFactor(TH2F *scaleMap, double pt, double eta, float sf_err, bool abs) {
//    std::cout<<"called getScaleFactor"<<std::endl;
  //  std::cout<<pt<<":, "<<eta<<std::endl;
    float sfactor = 1.0;
    int binx = scaleMap->GetXaxis()->FindBin(pt);
	 int biny;
    if (abs==0) biny = scaleMap->GetYaxis()->FindBin(eta);
    else biny = scaleMap->GetYaxis()->FindBin(TMath::Abs(eta));
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

using namespace std;

void TMVAapply_all_syst::Loop(TString inputfile, TString output_dir, TString region)
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
	TFile *output = TFile::Open(output_dir+"_"+region+".root","recreate");
	TTree *tree = fChain->CloneTree(0);
	const int Nsyst=5;
	float BDT_VBF[Nsyst];
	int pass_sel[Nsyst];
//	TBranch *branchBDT_VBF = tree->Branch("BDT_VBF",&BDT_VBF,"BDT_VBF/F");
	TBranch *branchBDT_VBF[Nsyst];
	TString BDT_VBF_names[Nsyst] = {"BDT_VBF_JES_up","BDT_VBF_JES_down","BDT_VBF_JER_up","BDT_VBF_JER_down","BDT_VBF"};
	TBranch *branchPassSel[Nsyst];
	TBranch *branchWeightEff;
	float WeightEff = 1.;
	TBranch *branchWeightQGL;
	float WeightQGL = 1.;
//	TBranch *branchWeightDY_zpt;
	float WeightDY_zpt = 1.;
	TString pass_sel_names[Nsyst] = {"PassSelection_JES_up","PassSelection_JES_down","PassSelection_JER_up","PassSelection_JER_down","PassSelection_nom"};
	for (int je_type=0;je_type<Nsyst;je_type++){
		branchBDT_VBF[je_type] = tree->Branch(BDT_VBF_names[je_type],&BDT_VBF[je_type],BDT_VBF_names[je_type]+"/F");
		branchPassSel[je_type] = tree->Branch(pass_sel_names[je_type],&pass_sel[je_type],pass_sel_names[je_type]+"/I");
	}
		branchWeightEff = tree->Branch("WeightEff",&WeightEff,"WeightEff/F");
		branchWeightQGL = tree->Branch("WeightQGL",&WeightQGL,"WeightQGL/F");
//		branchWeightDY_zpt = tree->Branch("WeightDY_zpt",&WeightDY_zpt,"WeightDY_zpt/F");

	TString weightfile[2]= {"/mnt/t3nfs01/data01/shome/nchernya/VBFZll/main_mva/weights/TMVAClassification_BDTG_allmu_new_v25/TMVAClassification_BDTG.weights.xml", "/mnt/t3nfs01/data01/shome/nchernya/VBFZll/main_mva/weights/TMVAClassification_BDTG_allel_new_v25/TMVAClassification_BDTG.weights.xml"};
   TMVA::Reader *reader = new TMVA::Reader("Silent");
	float var1,var2,var3,var4,var5,var6,var7,var8,var9,var10, var11, var12,var13;
	reader->AddVariable("Mqq",&var1);
	reader->AddVariable("DeltaEtaQQ",&var2);
	reader->AddVariable("qgl_1q",&var3);
	reader->AddVariable("qgl_2q",&var4);
//	reader->AddVariable("Jet2q_pt",&var3);
//	reader->AddVariable("axis2_jet1",&var4);
//	reader->AddVariable("axis2_jet2",&var5);
	reader->AddVariable("qq_pt",&var5);
	reader->AddVariable("RptHard",&var6);
	reader->AddVariable("Zll_zstar",&var7);
//	reader->AddVariable("axis2_jet1",&var3);
//	reader->AddVariable("axis2_jet2",&var4);
//	reader->AddVariable("Zll_zstar",&var5);
	if (region.CompareTo("mu")==0) reader->BookMVA("BDTG", weightfile[0]);
	if (region.CompareTo("el")==0) reader->BookMVA("BDTG", weightfile[1]);

	TFile *input_file = TFile::Open(inputfile);
	TH1F*	Count = (TH1F*)input_file->Get("Count");
	TH1F*	CountPosWeight = (TH1F*)input_file->Get("CountPosWeight");
	TH1F*	CountNegWeight =(TH1F*)input_file->Get("CountNegWeight");
	TH1F*	CountWeighted =(TH1F*)input_file->Get("CountWeighted");
	TH1F* CountFullWeighted  = (TH1F*)input_file->Get("CountFullWeighted");
	TH1F* CountWeightedLHEWeightScale  = (TH1F*)input_file->Get("CountWeightedLHEWeightScale");
	TH1F* CountWeightedLHEWeightPdf  = (TH1F*)input_file->Get("CountWeightedLHEWeightPdf");
	TH1F* CountWeightedLHEWeightaTGC  = (TH1F*)input_file->Get("CountWeightedLHEWeightaTGC");
	TH1F* CountWeightedLHEWeightaTGC_nom  = (TH1F*)input_file->Get("CountWeightedLHEWeightaTGC_nom");
/*

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

	TFile* file_id_el = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/v25/TriggerEffMap_ScaleFactor_MVAIDWP80_80x.root");
	TH2F* id_el = (TH2F*)file_id_el->Get("TriggerEffMap_ScaleFactor_MVAIDWP80_80x");
	TFile* file_tracker_el = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/v25/TriggerEffMap_ScaleFactor_tracker_80x.root");
	TH2F* tracker_el = (TH2F*)file_tracker_el->Get("TriggerEffMap_ScaleFactor_tracker_80x");

*/


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


	
	int events_saved=0;

	float weight;	
 

	cout<<nentries<<endl;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
//   for (Long64_t jentry=0; jentry<200;jentry++) {
	   Long64_t ientry = LoadTree(jentry);
	   if (ientry < 0) break;
	   nb = fChain->GetEntry(jentry);   nbytes += nb;
			  
		if (isData==1) genWeight = 1.;
		if (json!=1) continue;

	

///////////////////////////SYStematics loop////////////////////
		bool passes_any=false;	
     for (int current_syst=0;current_syst<Nsyst;current_syst++){
		pass_sel[current_syst] = 0;	
		if ((isData==1) && (current_syst!=4)) continue;

		WeightEff=1.;
		WeightQGL=1.;
		WeightDY_zpt=1.;

//		Float_t Jet_pt_new[njets];
//		Float_t Jet_mass_new[njets];
		for (int i=0;i<nJet;i++){
			if (current_syst==0) {Jet_pt[i] = Jet_pt[i]*Jet_corr_JECUp[i]/Jet_corr[i];  Jet_mass[i]=Jet_mass[i]*Jet_corr_JECUp[i]/Jet_corr[i];  }
			if (current_syst==1) {Jet_pt[i] = Jet_pt[i]*Jet_corr_JECDown[i]/Jet_corr[i];Jet_mass[i]=Jet_mass[i]*Jet_corr_JECDown[i]/Jet_corr[i];}
			if (current_syst==2) {
				if ((Jet_corr_JERUp[i]<=0)||(Jet_corr_JER[i]<=0)){ Jet_pt[i] = Jet_pt[i]; Jet_mass[i] = Jet_mass[i];}
				else	{
					Jet_pt[i] =  Jet_pt[i]*Jet_corr_JERUp[i]/Jet_corr_JER[i]; 
					Jet_mass[i]=Jet_mass[i]*Jet_corr_JERUp[i]/Jet_corr_JER[i];
				}
			}
			if (current_syst==3) {
				if ((Jet_corr_JERDown[i]<=0)||(Jet_corr_JER[i]<=0)){ Jet_pt[i] = Jet_pt[i]; Jet_mass[i] = Jet_mass[i];}
				else {
					Jet_pt[i] =  Jet_pt[i]*Jet_corr_JERDown[i]/Jet_corr_JER[i];
					Jet_mass[i]=Jet_mass[i]*Jet_corr_JERDown[i]/Jet_corr_JER[i];
				} 
			}
		}


		
		double GenVbosons_pt_first = GenVbosons_pt[0];
		int GenVbosons_pdgId_first = GenVbosons_pdgId[0];


		int good_jets = 0;
		int pt_num1 = -1;
		int pt_num2 = -1;
		TLorentzVector Qjet1;
		TLorentzVector Qjet2;
		TLorentzVector qq;
		vector<TLorentzVector> jets_pv;
		vector<int> jets_indices;
		///////////////////////
		//preselection/////
		for (int i=0;i<nJet;i++){
			TLorentzVector jet0;
			if (!((Jet_id[i]>2)&&(Jet_puId[i]>0))) continue;
			jet0.SetPtEtaPhiM(Jet_pt[i],Jet_eta[i],Jet_phi[i],Jet_mass[i]);
			jets_pv.push_back(jet0);
			jets_indices.push_back(i);
			good_jets++;
		}
		if (good_jets<2) continue;

		Qjet1 = jets_pv[0];
		Qjet2 = jets_pv[1];
		float jet3_pt = 0. ;
		if (good_jets>=3) jet3_pt = jets_pv[2].Pt();
	//	else jet3_pt = 10.;
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
			if (!((selLeptons_eleMVAIdSppring16GenPurp[i]>=2)&& (selLeptons_relIso03[i]<0.15) && (TMath::Abs(selLeptons_pdgId[i])==11))) continue;
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
	//		if ((count_l==1) && (selLeptons_charge[idx_1stLepton]*selLeptons_charge[i] > 0)) continue;
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
		if ((selLeptons_charge[idx_1stLepton]*selLeptons_charge[idx_2ndLepton])>0) continue;

///////////////muon corrections 2016 calibration////////////////
		if (region.CompareTo("mu")==0) {
			double dataSF1, dataSF2 ;
			double mcSF1, mcSF2;
			if (isData==1) {
				dataSF1 = rc->kScaleDT(selLeptons_charge[idx_1stLepton], lepton1.Pt(), lepton1.Eta(), lepton1.Phi(), 0, 0);
				dataSF2 = rc->kScaleDT(selLeptons_charge[idx_2ndLepton], lepton2.Pt(), lepton2.Eta(), lepton2.Phi(), 0, 0);
				lepton1.SetPtEtaPhiM(lepton1.Pt()*dataSF1,lepton1.Eta(), lepton1.Phi(), lepton1.M() );
				lepton2.SetPtEtaPhiM(lepton2.Pt()*dataSF2,lepton2.Eta(), lepton2.Phi(), lepton2.M() );
			}
			if (isData!=1) {
				double u1 = gRandom->Rndm();
				double u2 = gRandom->Rndm();
				mcSF1 = rc->kScaleAndSmearMC(selLeptons_charge[idx_1stLepton], lepton1.Pt(), lepton1.Eta(), lepton1.Phi(),  selLeptons_trackerLayers[idx_1stLepton], u1, u2, 0, 0);
				mcSF2 = rc->kScaleAndSmearMC(selLeptons_charge[idx_2ndLepton], lepton2.Pt(), lepton2.Eta(), lepton2.Phi(),  selLeptons_trackerLayers[idx_2ndLepton], u1, u2, 0, 0);
				lepton1.SetPtEtaPhiM(lepton1.Pt()*mcSF1,lepton1.Eta(), lepton1.Phi(), lepton1.M() );
			}
		}

////////////////////////////////////////////////////////////////
			if  (region.CompareTo("mu")==0) if (!((HLT_BIT_HLT_IsoMu24_v==1) || (HLT_BIT_HLT_IsoTkMu24_v==1)  )) continue; 
			if  (region.CompareTo("el")==0) if (!(HLT_BIT_HLT_Ele27_eta2p1_WPTight_Gsf_v == 1)) continue;

////////////////////////////////////////////////////////////////
		if (isData!=1) {
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

				WeightEff= eff1*eff1_id*eff2_id*eff1_iso*eff2_iso*eff1_tracker*eff2_tracker; 	
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
				WeightEff*=eff1* eff1_id*eff2_id*eff1_tracker*eff2_tracker; 	
			}
		}

//		string file_tag_str = file_tag.Data();
	//	if  (file_tag_str.find("DYJetstoLL")!=std::string::npos)  genweight*=ptWeightEWK(nGenVbosons, GenVbosons_pt_first, VtypeSim, GenVbosons_pdgId_first); 

//////////////////////QGL///////////////////////////
WeightQGL=1.;
		if (!( (isData==1)|| (Jet_partonFlavour[jets_indices[0]] ==0 ) || (TMath::Abs(Jet_eta[jets_indices[0]])>=2) || (Jet_qgl[jets_indices[0]] < 0) ) ) {
			if (TMath::Abs(Jet_partonFlavour[jets_indices[0]]) < 4 ) WeightQGL*=func_qgl_q->Eval(Jet_qgl[jets_indices[0]]);
			if (TMath::Abs(Jet_partonFlavour[jets_indices[0]]) ==21 ) WeightQGL*=func_qgl_g->Eval(Jet_qgl[jets_indices[0]]);
		}
	//	cout<<WeightQGL<<endl;
		if (!( (isData==1)|| (Jet_partonFlavour[jets_indices[1]] ==0 ) || (TMath::Abs(Jet_eta[jets_indices[1]])>=2) || (Jet_qgl[jets_indices[1]] < 0)) ) {
			if (TMath::Abs(Jet_partonFlavour[jets_indices[1]]) < 4 ) WeightQGL*=func_qgl_q->Eval(Jet_qgl[jets_indices[1]]);
			if (TMath::Abs(Jet_partonFlavour[jets_indices[1]]) ==21 ) WeightQGL*=func_qgl_g->Eval(Jet_qgl[jets_indices[1]]);
		}
			


//////////////////////////////////////////////////



		Zll = lepton1+lepton2;	
		Float_t Zll_mass = Zll.M();
		Float_t Zll_pt = Zll.Pt();
		Float_t Zll_eta = Zll.Eta();
		Float_t Zll_phi = Zll.Phi();

		if (good_jets<2) continue;
		if (Qjet1.Pt() < 50) continue;
		if (Qjet2.Pt() < 30) continue;
		if (Mqq<200) continue;
	
		if (lepton1.Pt()<30) continue;	
		if (lepton2.Pt()<20) continue;	
	 	if (region.CompareTo("el")==0) {
			if (TMath::Abs(lepton1.Eta())>2.1) continue;	
			if (TMath::Abs(lepton2.Eta())>2.1) continue;
		}	
	 	if (region.CompareTo("mu")==0) {
			if (TMath::Abs(lepton1.Eta())>2.4) continue;	
			if (TMath::Abs(lepton2.Eta())>2.4) continue;
		}	
		if (Zll_mass < 50 ) continue;
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
///////////////PUT it BACK, that was for VBFHmumu////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
		if (TMath::Abs(Zll_mass - 91.1876)>15) continue;

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
		passes_any=true;
		pass_sel[current_syst]=1;

		float Zll_ystar = Zll.Rapidity() - (Qjet1.Rapidity() + Qjet2.Rapidity()) ;
		float zll_zstar = TMath::Abs( Zll_ystar/ (Qjet1.Rapidity()-Qjet2.Rapidity() )) ;
		Float_t DeltaEtaQQSum = TMath::Abs(Qjet1.Eta()) +  TMath::Abs(Qjet2.Eta());
		Float_t PhiZQ1 = TMath::Abs(Zll.DeltaPhi(Qjet1));
		Float_t DeltaRelQQ = (Qjet1+Qjet2).Pt()/( Qjet1.Pt()+Qjet2.Pt()); 
		Float_t RptHard = (Qjet1+Qjet2+ Zll).Pt()/( Qjet1.Pt()+Qjet2.Pt() + Zll.Pt()); 




		var1 = Mqq;
		var2 = qqDeltaEta; 
//		var3 = Qjet2.Pt(); 
//		var4=TMath::Exp((-1)*Jet_axis2[jets_indices[0]]);
//		var5=TMath::Exp((-1)*Jet_axis2[jets_indices[1]]);
		var3=Jet_qgl[jets_indices[0]];
		var4=Jet_qgl[jets_indices[1]];
		var5 = qq_pt;
		var6= RptHard;
		var7 = zll_zstar;
//		var3=Jet_qgl[jets_indices[0]];
//		var4=Jet_qgl[jets_indices[1]];
//		var3=TMath::Exp((-1)*Jet_axis2[jets_indices[0]]);
//		var4=TMath::Exp((-1)*Jet_axis2[jets_indices[1]]);
//		var5 = zll_zstar;
	//	var5 = qq_pt;
//		var6= RptHard;
		BDT_VBF[current_syst] = reader->EvaluateMVA("BDTG");
	  } //// systematics loop ends
		if (passes_any==true)	tree->Fill();
//		for (int i=0;i<Nsyst;i++)
//			cout<<pass_sel[i]<<"  ";
//		cout<<endl;


	}  
	delete reader;
	output->cd();
	tree->AutoSave();
	Count->Write();
	CountPosWeight->Write();
	CountNegWeight->Write();
	CountWeighted->Write();
	CountFullWeighted->Write();
	CountWeightedLHEWeightScale->Write();
	CountWeightedLHEWeightPdf->Write();
	CountWeightedLHEWeightaTGC->Write();
	CountWeightedLHEWeightaTGC_nom->Write();
	output->Close();
	input_file->Close();


}

