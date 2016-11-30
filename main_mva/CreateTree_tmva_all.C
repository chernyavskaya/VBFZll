#define CreateTree_tmva_all_cxx
#include "CreateTree_tmva_all.h"
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

void CreateTree_tmva_all::Loop(TString file_tag, TString region)
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
	TMVAstruct TMVA;

gROOT->ProcessLine(".L /mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/RoccoR.cc++");
gROOT->ProcessLine(".L /mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/rochcor2016.cc++");

	TF1* func_JetsPt = new TF1("func_JetsPt","pol5",4.3,10);

	TF1* func_Mqq = new TF1("func_Mqq","pol6",0,10);
	func_Mqq->FixParameter(0,5968.223851);
	func_Mqq->FixParameter(1,-5554.340558);
	func_Mqq->FixParameter(2,2146.273308);
	func_Mqq->FixParameter(3,-440.733829);
	func_Mqq->FixParameter(4,50.729466);
	func_Mqq->FixParameter(5,-3.103351);
	func_Mqq->FixParameter(6,0.0788278);

	rochcor2016 *rmcor = new rochcor2016();
	TFile* file_trig_el = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/skimmed/TriggerEffMap_electronTriggerEfficiencyHLT_Ele27_WPLoose_eta2p1_WP90_BCDEF.root");
	TH2F* trig_el = (TH2F*)file_trig_el->Get("TriggerEffMap_electronTriggerEfficiencyHLT_Ele27_WPLoose_eta2p1_WP90_BCDEF");
	TFile* file_trig_mu_bf = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/skimmed/TriggerEffMap_MuonTrigger_data_all_IsoMu22_OR_IsoTkMu22_pteta_Run2016B_beforeL2Fix.root");
	TH2F* trig_mu_bf = (TH2F*)file_trig_mu_bf->Get("TriggerEffMap_MuonTrigger_data_all_IsoMu22_OR_IsoTkMu22_pteta_Run2016B_beforeL2Fix");
	TFile* file_trig_mu_aft = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/skimmed/TriggerEffMap_MuonTrigger_data_all_IsoMu22_OR_IsoTkMu22_pteta_Run2016B_afterL2Fix.root");
	TH2F* trig_mu_aft = (TH2F*)file_trig_mu_aft->Get("TriggerEffMap_MuonTrigger_data_all_IsoMu22_OR_IsoTkMu22_pteta_Run2016B_afterL2Fix");

	int events_saved=0;

	float weight;	
 

	TFile file("main_tmva_tree_"+file_tag+"_v24"+region+"_MqqLog.root","recreate");
	TTree *tree0 = new TTree("TMVA","TMVA");

	tree0->Branch("Mqq",&TMVA.Mqq,"Mqq/F");
	tree0->Branch("DeltaEtaQQ",&TMVA.DeltaEtaQQ,"DeltaEtaQQ/F");
	tree0->Branch("Jet2q_pt",&TMVA.Jet2q_pt,"Jet2q_pt/F");
	tree0->Branch("Jet1q_pt",&TMVA.Jet1q_pt,"Jet1q_pt/F");
	tree0->Branch("Jet1q_leadTrackPt",&TMVA.Jet1q_leadTrackPt,"Jet1q_leadTrackPt/F");
	tree0->Branch("Jet2q_leadTrackPt",&TMVA.Jet2q_leadTrackPt,"Jet2q_leadTrackPt/F");
	tree0->Branch("axis2_jet1",&TMVA.axis2_jet1,"axis2_jet1/F");
	tree0->Branch("axis2_jet2",&TMVA.axis2_jet2,"axis2_jet2/F");
	tree0->Branch("qq_pt",&TMVA.qq_pt,"qq_pt/F");
	tree0->Branch("Jet3_pt",&TMVA.Jet3_pt,"Jet3_pt/F");
	tree0->Branch("RptHard",&TMVA.RptHard,"RptHard/F");
	tree0->Branch("Zll_pt",&TMVA.Zll_pt,"Zll_pt/F");
	tree0->Branch("Zll_zstar",&TMVA.Zll_zstar,"Zll_zstar/F");
	tree0->Branch("weight",&weight);


	
	cout<<nentries<<endl;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
	   Long64_t ientry = LoadTree(jentry);
	   if (ientry < 0) break;
	   nb = fChain->GetEntry(jentry);   nbytes += nb;
			  
		if (isData==1) genWeight = 1.;
		if (genWeight <0) continue;
		if (json!=1) continue;

		
		if (region.CompareTo("mu")==0) if (!(Vtype==0)) continue;
		if (region.CompareTo("el")==0) if (!(Vtype==1)) continue;
		
		double GenVbosons_pt_first = GenVbosons_pt[0];
		int GenVbosons_pdgId_first = GenVbosons_pdgId[0];
		if  ((file_tag.CompareTo("DYJetstoLL")==0) 
|| (file_tag.CompareTo("DYJetstoLL_amc")==0) || (file_tag.CompareTo("DYJetstoLL_Pt-100_amc")==0)|| (file_tag.CompareTo("DYJetstoLL_Pt-100To250_amc")==0)|| (file_tag.CompareTo("DYJetstoLL_Pt-250To400_amc")==0) || (file_tag.CompareTo("DYJetstoLL_Pt-400To650_amc")==0)|| (file_tag.CompareTo("DYJetstoLL_Pt-650ToInf_amc")==0)  
 || (file_tag.CompareTo("DYJetstoLL_HT100")==0) ||(file_tag.CompareTo("DYJetstoLL_HT100_200")==0) ||(file_tag.CompareTo("DYJetstoLL_HT200_400")==0) || (file_tag.CompareTo("DYJetstoLL_HT400_600")==0) ||(file_tag.CompareTo("DYJetstoLL_HT600_Inf")==0)   ) genWeight*=ptWeightEWK(nGenVbosons, GenVbosons_pt_first, VtypeSim, GenVbosons_pdgId_first); 


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
			string file_tag_str = file_tag.Data();
			if  (file_tag_str.find("SingleMuon")!=std::string::npos) if (!((HLT_BIT_HLT_IsoMu27_v==1) || (HLT_BIT_HLT_IsoTkMu27_v==1) )) continue; 
			if  (file_tag.CompareTo("SingleElectron")==0) if (!(HLT_BIT_HLT_Ele27_eta2p1_WPLoose_Gsf_v == 1)) continue;
		} else if (isData!=1) {
			if (region.CompareTo("mu")==0) {
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
			}			
		}
		Float_t Mqq_log = TMath::Log(Mqq);	
		float jets_ptSum =TMath::Log(jets_pv[0].Pt() + jets_pv[1].Pt());
		if (region.CompareTo("el")==0) {
			func_JetsPt->FixParameter(0,-467.774);
			func_JetsPt->FixParameter(1,408.471);
			func_JetsPt->FixParameter(2,-141.63);
			func_JetsPt->FixParameter(3,24.4376);
			func_JetsPt->FixParameter(4,-2.0985);
			func_JetsPt->FixParameter(5,0.0717106);


			if ((Mqq_log<8.176 )&&(Mqq_log>5.256)) if (file_tag.CompareTo("DYJetstoLL")==0) genWeight*=func_Mqq->Eval(Mqq_log);		
//			if ((jets_ptSum>=4.3) && (jets_ptSum<=7.7113)) if (file_tag.CompareTo("DYJetstoLL")==0)  genWeight*=func_JetsPt->Eval(TMath::Log(jets_pv[0].Pt() + jets_pv[1].Pt()));		
		}
		if (region.CompareTo("mu")==0) {
			func_JetsPt->FixParameter(0,-328.522);
			func_JetsPt->FixParameter(1,280.199 );
			func_JetsPt->FixParameter(2,-94.6805);
			func_JetsPt->FixParameter(3,15.9021);
			func_JetsPt->FixParameter(4,-1.32777);
			func_JetsPt->FixParameter(5,0.0440635);

	//		if ((jets_ptSum>=4.3) && (jets_ptSum<=8.145)) if (file_tag.CompareTo("DYJetstoLL")==0)   genWeight*=func_JetsPt->Eval(TMath::Log(jets_pv[0].Pt() + jets_pv[1].Pt()));		
			if ((Mqq_log<8.176 )&&(Mqq_log>5.256)) if (file_tag.CompareTo("DYJetstoLL")==0)  genWeight*=func_Mqq->Eval(Mqq_log);		
		}


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




		TMVA.Mqq = Mqq;
		TMVA.Jet2q_pt = Qjet2.Pt();
		if (zll_zstar<5) TMVA.Zll_zstar = zll_zstar;
		TMVA.axis2_jet1 = Jet_axis2[jets_indices[0]];
		TMVA.Jet1q_pt = Qjet1.Pt();
		TMVA.axis2_jet2 = Jet_axis2[jets_indices[1]];
		TMVA.Jet1q_leadTrackPt = Jet_leadTrackPt[jets_indices[0]];
		TMVA.RptHard = RptHard;
		TMVA.DeltaEtaQQ = qqDeltaEta;
		TMVA.Jet2q_leadTrackPt = Jet_leadTrackPt[jets_indices[1]];
		TMVA.qq_pt = qq_pt;
		TMVA.Zll_pt = Zll_pt;
		TMVA.Jet3_pt = jet3_pt;
		weight=genWeight/TMath::Abs(genWeight)*puWeight;
		

		tree0->Fill();	
		events_saved++;		
		if (region.CompareTo("mu")==0) if (events_saved>=140000) break;
		if (region.CompareTo("el")==0) if (events_saved>=80000) break;

	}  

	file.Write();
	file.Close();

}

