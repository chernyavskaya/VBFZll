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


using namespace std;

void CreateTree_tmva_all::Loop(TString file_tag, TString region)
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
	TMVAstruct TMVA;


	int events_saved=0;

	float weight;	
 

	TFile file("main_tmva_tree_"+file_tag+"_v24"+region+".root","recreate");
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
		else jet3_pt = 10.;
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
			if  (file_tag_str.find("SingleMuon")!=std::string::npos) if (!((HLT_BIT_HLT_IsoMu22_v==1) || (HLT_BIT_HLT_IsoTkMu22_v==1) )) continue; 
			if  (file_tag.CompareTo("SingleElectron")==0) if (!(HLT_BIT_HLT_Ele27_eta2p1_WPLoose_Gsf_v == 1)) continue;
		} else if (isData!=1) {
			if (region.CompareTo("mu")==0) {
			//	cout<<vLeptons_SF_HLT_RunD4p2[0]<<"  "<<vLeptons_SF_HLT_RunD4p3[0]<<endl;
				genWeight*=(0.032*vLeptons_SF_HLT_RunD4p2[0] + 0.96799*vLeptons_SF_HLT_RunD4p3[0]);
			//	cout<<genWeight<<endl;
				genWeight*= vLeptons_SF_IdCutLoose[0]*vLeptons_SF_IdCutLoose[1] * vLeptons_SF_IsoLoose[0]* vLeptons_SF_IsoLoose[1]* vLeptons_SF_trk_eta[0]*vLeptons_SF_trk_eta[1];
			}
		}

		if (good_jets<2) continue;
		if (Qjet1.Pt() < 50) continue;
		if (Qjet1.Pt() < 50) continue;
		if (Qjet2.Pt() < 30) continue;
		if (Mqq<200) continue;
	
		if (lepton1.Pt()<20) continue;	
		if (lepton2.Pt()<20) continue;	
	 	if (region.CompareTo("el")==0) {
			if (TMath::Abs(lepton1.Eta())>2.1) continue;	
			if (TMath::Abs(lepton2.Eta())>2.1) continue;
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

	}  

	file.Write();
	file.Close();

}

