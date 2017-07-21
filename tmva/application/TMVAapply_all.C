#define TMVAapply_all_cxx
#include "TMVAapply_all.h"
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
#include "/mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/rochcor2016.h"
#include "/mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/rochcor2016.cc"
#include "/mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/RoccoR.cc"
#include "/mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/RoccoR.h"


const int njets = 300;

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

void TMVAapply_all::Loop(TString inputfile, TString output_dir, TString region)
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
	TFile *output = TFile::Open(output_dir+"_"+region+".root","recreate");
	TTree *tree = fChain->CloneTree(0);
	float BDT_VBF;
	TBranch *branchBDT_VBF = tree->Branch("BDT_VBF",&BDT_VBF,"BDT_VBF/F");
	TString weightfile[2]= {"/mnt/t3nfs01/data01/shome/nchernya/VBFZll/main_mva/weights/TMVAClassification_BDTG_allmu_v24/TMVAClassification_BDTG.weights.xml", "/mnt/t3nfs01/data01/shome/nchernya/VBFZll/main_mva/weights/TMVAClassification_BDTG_allel_v24/TMVAClassification_BDTG.weights.xml"};
   TMVA::Reader *reader = new TMVA::Reader("Silent");
	float var1,var2,var3,var4,var5,var6,var7,var8,var9,var10, var11, var12,var13;
	reader->AddVariable("Mqq",&var1);
	reader->AddVariable("DeltaEtaQQ",&var2);
	reader->AddVariable("Jet2q_pt",&var3);
	reader->AddVariable("axis2_jet1",&var4);
	reader->AddVariable("axis2_jet2",&var5);
	reader->AddVariable("qq_pt",&var6);
	reader->AddVariable("RptHard",&var7);
	reader->AddVariable("Zll_zstar",&var8);
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

gROOT->ProcessLine(".L /mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/RoccoR.cc++");
gROOT->ProcessLine(".L /mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/rochcor2016.cc++");

//	TF1* func_JetsPt = new TF1("func_JetsPt","pol5",4.3,10);
	rochcor2016 *rmcor = new rochcor2016();
//	TFile* file_trig_el = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/skimmed/TriggerEffMap_electronTriggerEfficiencyHLT_Ele27_WPLoose_eta2p1_WP90_BCDEF.root");
//	TH2F* trig_el = (TH2F*)file_trig_el->Get("TriggerEffMap_electronTriggerEfficiencyHLT_Ele27_WPLoose_eta2p1_WP90_BCDEF");
//	TFile* file_trig_mu_bf = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/skimmed/TriggerEffMap_MuonTrigger_data_all_IsoMu22_OR_IsoTkMu22_pteta_Run2016B_beforeL2Fix.root");
//	TH2F* trig_mu_bf = (TH2F*)file_trig_mu_bf->Get("TriggerEffMap_MuonTrigger_data_all_IsoMu22_OR_IsoTkMu22_pteta_Run2016B_beforeL2Fix");
//	TFile* file_trig_mu_aft = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/skimmed/TriggerEffMap_MuonTrigger_data_all_IsoMu22_OR_IsoTkMu22_pteta_Run2016B_afterL2Fix.root");
//	TH2F* trig_mu_aft = (TH2F*)file_trig_mu_aft->Get("TriggerEffMap_MuonTrigger_data_all_IsoMu22_OR_IsoTkMu22_pteta_Run2016B_afterL2Fix");

	int events_saved=0;

	float weight;	
 

	cout<<nentries<<endl;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
	   Long64_t ientry = LoadTree(jentry);
	   if (ientry < 0) break;
	   nb = fChain->GetEntry(jentry);   nbytes += nb;
			  
		if (isData==1) genWeight = 1.;
		if (json!=1) continue;

		
		if (region.CompareTo("mu")==0) if (!(Vtype==0)) continue;
		if (region.CompareTo("el")==0) if (!(Vtype==1)) continue;
		
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
		lepton1.SetPtEtaPhiM(vLeptons_pt[0], vLeptons_eta[0], vLeptons_phi[0], vLeptons_mass[0]);	
		int idx_2ndLepton = 0;
		for (int i=1; i<nvLeptons;i++ ){
			if (vLeptons_charge[0]*vLeptons_charge[i] < 0) {
				idx_2ndLepton=i;
				break;
			}
		}
		lepton2.SetPtEtaPhiM(vLeptons_pt[idx_2ndLepton], vLeptons_eta[idx_2ndLepton], vLeptons_phi[idx_2ndLepton], vLeptons_mass[idx_2ndLepton]);
		Zll = lepton1+lepton2;	
		Float_t Zll_mass = Zll.M();
		Float_t Zll_pt = Zll.Pt();
		Float_t Zll_eta = Zll.Eta();
		Float_t Zll_phi = Zll.Phi();


		if (isData==1) { 
			if  (region.CompareTo("mu")==0) if (!((HLT_BIT_HLT_IsoMu27_v==1) || (HLT_BIT_HLT_IsoTkMu27_v==1) )) continue; 
			if  (region.CompareTo("el")==0) if (!(HLT_BIT_HLT_Ele27_eta2p1_WPLoose_Gsf_v == 1)) continue;
		} else if (isData!=1) {
		/*	if (region.CompareTo("mu")==0) {
				float SF_mu_bf_err1 = 0.;
				float SF_mu_bf_err2 = 0.;
				float SF_mu_aft_err1 = 0.;
				float SF_mu_aft_err2 = 0.;
				bool abs=1;
				float eff1 =0.02772*getScaleFactor(trig_mu_bf, lepton1.Pt(), lepton1.Eta(), SF_mu_bf_err1,abs ) + 0.97227*getScaleFactor(trig_mu_aft, lepton1.Pt(), lepton1.Eta(), SF_mu_bf_err1,abs ) ;  	
				float eff2 =0.02772*getScaleFactor(trig_mu_bf, lepton2.Pt(), lepton2.Eta(), SF_mu_bf_err2,abs ) + 0.97227*getScaleFactor(trig_mu_aft, lepton2.Pt(), lepton2.Eta(), SF_mu_bf_err2,abs  ) ;  
				genWeight*= eff1*(1-eff2)*eff1 + eff2*(1-eff1)*eff2 + eff1*eff1*eff2*eff2; 	
				genWeight*= vLeptons_SF_IdCutLoose[0]*vLeptons_SF_IdCutLoose[1] * vLeptons_SF_IsoLoose[0]* vLeptons_SF_IsoLoose[1]* vLeptons_SF_trk_eta[0]*vLeptons_SF_trk_eta[1];
			}
			if (region.CompareTo("el")==0) {
				float SF_el_err1 = 0.;
				float SF_el_err2 = 0.;
				bool abs=0;
				float eff1 = getScaleFactor(trig_el, lepton1.Pt(), lepton1.Eta(), SF_el_err1,abs );  	
				float eff2 = getScaleFactor(trig_el, lepton2.Pt(), lepton2.Eta(), SF_el_err2,abs ); 
				genWeight*= eff1*(1-eff2)*eff1 + eff2*(1-eff1)*eff2 + eff1*eff1*eff2*eff2; 	
			}	*/		
		}
/*		float jets_ptSum =TMath::Log(jets_pv[0].Pt() + jets_pv[1].Pt());
		if (region.CompareTo("el")==0) {
			func_JetsPt->FixParameter(0,-467.774);
			func_JetsPt->FixParameter(1,408.471);
			func_JetsPt->FixParameter(2,-141.63);
			func_JetsPt->FixParameter(3,24.4376);
			func_JetsPt->FixParameter(4,-2.0985);
			func_JetsPt->FixParameter(5,0.0717106);


			if ((jets_ptSum>=4.3) && (jets_ptSum<=7.7113)) if (file_tag.CompareTo("DYJetstoLL")==0)  genWeight*=func_JetsPt->Eval(TMath::Log(jets_pv[0].Pt() + jets_pv[1].Pt()));		
		}
		if (region.CompareTo("mu")==0) {
			func_JetsPt->FixParameter(0,-328.522);
			func_JetsPt->FixParameter(1,280.199 );
			func_JetsPt->FixParameter(2,-94.6805);
			func_JetsPt->FixParameter(3,15.9021);
			func_JetsPt->FixParameter(4,-1.32777);
			func_JetsPt->FixParameter(5,0.0440635);

			if ((jets_ptSum>=4.3) && (jets_ptSum<=8.145)) if (file_tag.CompareTo("DYJetstoLL")==0)   genWeight*=func_JetsPt->Eval(TMath::Log(jets_pv[0].Pt() + jets_pv[1].Pt()));		
		}
*/

		float qter1 = 1.0;
		float qter2 = 1.0;
		float mu_correction1 = 1.0;
		float mu_correction2 = 1.0;
		if (isData!=1) 	if (region.CompareTo("mu")==0) {
			rmcor->momcor_mc(lepton1, vLeptons_charge[0], vLeptons_trackerLayers[0], qter1);
			rmcor->momcor_mc(lepton2, vLeptons_charge[idx_2ndLepton], vLeptons_trackerLayers[idx_2ndLepton], qter2);
			}


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
		if (TMath::Abs(Zll_mass - 91.1876)>15) continue;

		float Zll_ystar = Zll.Rapidity() - (Qjet1.Rapidity() + Qjet2.Rapidity()) ;
		float zll_zstar = TMath::Abs( Zll_ystar/ (Qjet1.Rapidity()-Qjet2.Rapidity() )) ;
		Float_t DeltaEtaQQSum = TMath::Abs(Qjet1.Eta()) +  TMath::Abs(Qjet2.Eta());
		Float_t PhiZQ1 = TMath::Abs(Zll.DeltaPhi(Qjet1));
		Float_t DeltaRelQQ = (Qjet1+Qjet2).Pt()/( Qjet1.Pt()+Qjet2.Pt()); 
		Float_t RptHard = (Qjet1+Qjet2+ Zll).Pt()/( Qjet1.Pt()+Qjet2.Pt() + Zll.Pt()); 




		var1 = Mqq;
		var2 = qqDeltaEta; 
		var3 = Qjet2.Pt();
		var4=Jet_axis2[jets_indices[0]];
		var5=Jet_axis2[jets_indices[1]];
		var6= qq_pt;
		var7= RptHard;
		var8 = zll_zstar;
		BDT_VBF = reader->EvaluateMVA("BDTG");
		tree->Fill();


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
	output->Close();


}

