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
#include "/mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/mucorr/2016/RoccoR.cc"


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
	Float_t qgl_1q;
	Float_t qgl_2q;
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

void CreateTree_tmva_all::Loop(TString file_tag, TString region, TString file_name)
{
   if (fChain == 0) return;


std::map <TString, float> xsec;
xsec["EWK_LL_JJ"]=1.664;
xsec["DYJetstoLL_madgraph"] =  5765.4;
xsec["DYJetstoLL_HT100"] =  5765.4;
xsec["DYJetstoLL_HT100to200"] = 181.302; 
xsec["DYJetstoLL_HT200to400"] =50.4177  ; 
xsec["DYJetstoLL_HT400to600"] =6.98394; 
xsec["DYJetstoLL_HT600to800"] = 1.6814;
xsec["DYJetstoLL_HT800to1200"] = 0.7754;
xsec["DYJetstoLL_HT1200to2500"] = 0.186;
xsec["DYJetstoLL_HT2500toInf"] = 0.00438495;
 

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
	TMVAstruct TMVA;

	TFile *file_initial = TFile::Open(file_name);
 	TH1F *countWeighted = (TH1F*)file_initial->Get("CountWeighted");
 	TH1F *countPosWeighted = (TH1F*)file_initial->Get("CountPosWeight");
 	TH1F *countNegWeighted = (TH1F*)file_initial->Get("CountNegWeight");
	Int_t events_generated = countWeighted->GetBinContent(1);
	if (events_generated==0) events_generated =  countPosWeighted->GetBinContent(1) - countNegWeighted->GetBinContent(1);

	cout<< events_generated<<endl;

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
	

	RoccoR  *rc = new RoccoR("/mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/mucorr/2016/rcdata.2016.v3/");

	TF1* func_JetsPt = new TF1("func_JetsPt","pol5",4.3,10);

	TF1* func_Mqq = new TF1("func_Mqq","pol6",0,10);
	func_Mqq->FixParameter(0,-1378.361758);
	func_Mqq->FixParameter(1,1327.383178);
	func_Mqq->FixParameter(2,-529.261759);
	func_Mqq->FixParameter(3,111.854413);
	func_Mqq->FixParameter(4,-13.208941);
	func_Mqq->FixParameter(5,0.826140);
	func_Mqq->FixParameter(6,-0.021376);


	int events_saved=0;

	float weight;	
 

	TFile file("main_tmva_tree_"+file_tag+"_v25"+region+"_new.root","recreate");
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
	tree0->Branch("qgl_1q",&TMVA.qgl_1q,"qgl_1q/F");
	tree0->Branch("qgl_2q",&TMVA.qgl_2q,"qgl_2q/F");
	tree0->Branch("Jet3_pt",&TMVA.Jet3_pt,"Jet3_pt/F");
	tree0->Branch("RptHard",&TMVA.RptHard,"RptHard/F");
	tree0->Branch("Zll_pt",&TMVA.Zll_pt,"Zll_pt/F");
	tree0->Branch("Zll_zstar",&TMVA.Zll_zstar,"Zll_zstar/F");
	tree0->Branch("weight",&weight);


	
	cout<<nentries<<endl;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
 //  for (Long64_t jentry=0; jentry<2000;jentry++) {
	   Long64_t ientry = LoadTree(jentry);
	   if (ientry < 0) break;
	   nb = fChain->GetEntry(jentry);   nbytes += nb;
			  
		if (isData==1) genWeight = 1.;
		if (genWeight <0) continue;
		if (json!=1) continue;
		weight=1;

		weight=1./events_generated*xsec[file_tag]*1000.; 
		
	//	if (region.CompareTo("mu")==0) if (!(Vtype==0)) continue;
	//	if (region.CompareTo("el")==0) if (!(Vtype==1)) continue;
		
		double GenVbosons_pt_first = GenVbosons_pt[0];
		int GenVbosons_pdgId_first = GenVbosons_pdgId[0];

		string file_tag_str = file_tag.Data();
		if  ((file_tag.CompareTo("DYJetstoLL_madgraph")==0) && (lheHT > 100)) continue;		


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
			if (!((selLeptons_eleMVAIdSppring16GenPurp[i]>=2)&& (selLeptons_relIso03[i]<0.15)&& (TMath::Abs(selLeptons_pdgId[i])==11))) continue;
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

	



		if (region.CompareTo("mu")==0) {
			if (isData!=1) {
				double mcSF1, mcSF2;
				double u1 = gRandom->Rndm();
				double u2 = gRandom->Rndm();
				mcSF1 = rc->kScaleAndSmearMC(selLeptons_charge[idx_1stLepton], lepton1.Pt(), lepton1.Eta(), lepton1.Phi(),  selLeptons_trackerLayers[idx_1stLepton], u1, u2, 0, 0);
				mcSF2 = rc->kScaleAndSmearMC(selLeptons_charge[idx_2ndLepton], lepton2.Pt(), lepton2.Eta(), lepton2.Phi(),  selLeptons_trackerLayers[idx_2ndLepton], u1, u2, 0, 0);
				lepton1.SetPtEtaPhiM(lepton1.Pt()*mcSF1,lepton1.Eta(), lepton1.Phi(), lepton1.M() );
				lepton2.SetPtEtaPhiM(lepton2.Pt()*mcSF2,lepton2.Eta(), lepton2.Phi(), lepton2.M() );
			}
		}

	


		Zll = lepton1+lepton2;	
		Float_t Zll_mass = Zll.M();
		Float_t Zll_pt = Zll.Pt();
		Float_t Zll_eta = Zll.Eta();
		Float_t Zll_phi = Zll.Phi();

			if  (region.CompareTo("mu")==0) if (!((HLT_BIT_HLT_IsoMu24_v==1) || (HLT_BIT_HLT_IsoTkMu24_v==1)  )) continue; 
			if  (region.CompareTo("el")==0) if (!(HLT_BIT_HLT_Ele27_eta2p1_WPTight_Gsf_v == 1)) continue;

		if (isData!=1) {
			if (region.CompareTo("mu")==0) {
				float SF_mu_bf_err1 = 0.;
				float SF_mu_bf_err2 = 0.;
				float SF_mu_aft_err1 = 0.;
				float SF_mu_aft_err2 = 0.;
				bool abs=1;
				float eff1 =20.1/36.4*getScaleFactor(trig_mu_bf, lepton1.Pt(), lepton1.Eta(), SF_mu_bf_err1,abs ) + 16.3/36.4*getScaleFactor(trig_mu_aft, lepton1.Pt(), lepton1.Eta(), SF_mu_bf_err1,abs ) ;  	
				float eff2 =20.1/36.4*getScaleFactor(trig_mu_bf, lepton2.Pt(), lepton2.Eta(), SF_mu_bf_err2,abs ) + 16.3/36.4*getScaleFactor(trig_mu_aft, lepton2.Pt(), lepton2.Eta(), SF_mu_bf_err2,abs  ) ;  
				weight*= eff1*(1-eff2)*eff1 + eff2*(1-eff1)*eff2 + eff1*eff1*eff2*eff2; 	
				float eff1_id =20.1/36.4*getScaleFactor(id_mu_bf, lepton1.Pt(), lepton1.Eta(), SF_mu_bf_err1,abs ) + 16.3/36.4*getScaleFactor(id_mu_aft, lepton1.Pt(), lepton1.Eta(), SF_mu_bf_err1,abs ) ;  	
				float eff2_id =20.1/36.4*getScaleFactor(id_mu_bf, lepton2.Pt(), lepton2.Eta(), SF_mu_bf_err2,abs ) + 16.3/36.4*getScaleFactor(trig_mu_aft, lepton2.Pt(), lepton2.Eta(), SF_mu_bf_err2,abs  ) ;  
				float eff1_iso =20.1/36.4*getScaleFactor(iso_mu_bf, lepton1.Pt(), lepton1.Eta(), SF_mu_bf_err1,abs ) + 16.3/36.4*getScaleFactor(id_mu_aft, lepton1.Pt(), lepton1.Eta(), SF_mu_bf_err1,abs ) ;  	
				float eff2_iso =20.1/36.4*getScaleFactor(iso_mu_bf, lepton2.Pt(), lepton2.Eta(), SF_mu_bf_err2,abs ) + 16.3/36.4*getScaleFactor(trig_mu_aft, lepton2.Pt(), lepton2.Eta(), SF_mu_bf_err2,abs  ) ;  
				weight*= eff1_id*eff2_id*eff1_iso*eff2_iso; 	

			}
			if (region.CompareTo("el")==0) {
				float SF_el_err1 = 0.;
				float SF_el_err2 = 0.;
				bool abs=0;
				float eff1_id =getScaleFactor(id_el, lepton1.Pt(), lepton1.Eta(), SF_el_err1,abs ) ;  	
				float eff2_id =getScaleFactor(id_el, lepton2.Pt(), lepton2.Eta(), SF_el_err2,abs  ) ;  
				float eff1_tracker =getScaleFactor(tracker_el, lepton1.Pt(), lepton1.Eta(), SF_el_err1,abs ) ;  	
				float eff2_tracker =getScaleFactor(tracker_el, lepton2.Pt(), lepton2.Eta(), SF_el_err2,abs  ) ;  
				weight*= eff1_id*eff2_id*eff1_tracker*eff2_tracker; 	
			}			
		}
		if  (file_tag_str.find("DYJetstoLL")!=std::string::npos)  weight*=ptWeightEWK(nGenVbosons, GenVbosons_pt_first, VtypeSim, GenVbosons_pdgId_first);
 
		Float_t Mqq_log = TMath::Log(Mqq);	

	if  (file_tag_str.find("DYJetstoLL")!=std::string::npos)   weight*=func_Mqq->Eval(Mqq_log);		
	/*
		Float_t Mqq_log = TMath::Log(Mqq);	
		float jets_ptSum =TMath::Log(jets_pv[0].Pt() + jets_pv[1].Pt());
		if (region.CompareTo("el")==0) {
			func_JetsPt->FixParameter(0,-467.774);
			func_JetsPt->FixParameter(1,408.471);
			func_JetsPt->FixParameter(2,-141.63);
			func_JetsPt->FixParameter(3,24.4376);
			func_JetsPt->FixParameter(4,-2.0985);
			func_JetsPt->FixParameter(5,0.0717106);


			if ((Mqq_log<8.176 )&&(Mqq_log>5.256)) if (file_tag.CompareTo("DYJetstoLL")==0) weight*=func_Mqq->Eval(Mqq_log);		
//			if ((jets_ptSum>=4.3) && (jets_ptSum<=7.7113)) if (file_tag.CompareTo("DYJetstoLL")==0)  weight*=func_JetsPt->Eval(TMath::Log(jets_pv[0].Pt() + jets_pv[1].Pt()));		
		}
		if (region.CompareTo("mu")==0) {
			func_JetsPt->FixParameter(0,-328.522);
			func_JetsPt->FixParameter(1,280.199 );
			func_JetsPt->FixParameter(2,-94.6805);
			func_JetsPt->FixParameter(3,15.9021);
			func_JetsPt->FixParameter(4,-1.32777);
			func_JetsPt->FixParameter(5,0.0440635);

	//		if ((jets_ptSum>=4.3) && (jets_ptSum<=8.145)) if (file_tag.CompareTo("DYJetstoLL")==0)   weight*=func_JetsPt->Eval(TMath::Log(jets_pv[0].Pt() + jets_pv[1].Pt()));		
			if ((Mqq_log<8.176 )&&(Mqq_log>5.256)) if (file_tag.CompareTo("DYJetstoLL")==0)  weight*=func_Mqq->Eval(Mqq_log);		
		}

*/


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
		TMVA.Jet1q_pt = Qjet1.Pt();
		if (zll_zstar<5) TMVA.Zll_zstar = zll_zstar;
		if (TMath::Exp((-1)*Jet_axis2[jets_indices[0]]) < 0.2)  TMVA.axis2_jet1 = TMath::Exp((-1)*Jet_axis2[jets_indices[0]]);
		if (TMath::Exp((-1)*Jet_axis2[jets_indices[1]]) < 0.2) TMVA.axis2_jet2 = TMath::Exp((-1)*Jet_axis2[jets_indices[1]]);
		TMVA.Jet1q_leadTrackPt = Jet_leadTrackPt[jets_indices[0]];
		TMVA.RptHard = RptHard;
		TMVA.DeltaEtaQQ = qqDeltaEta;
		TMVA.Jet2q_leadTrackPt = Jet_leadTrackPt[jets_indices[1]];
		TMVA.qq_pt = qq_pt;
		TMVA.Zll_pt = Zll_pt;
		TMVA.Jet3_pt = jet3_pt;
		TMVA.qgl_1q = Jet_qgl[jets_indices[0]];
		TMVA.qgl_2q = Jet_qgl[jets_indices[1]];
		weight=weight*genWeight/TMath::Abs(genWeight)*puWeight;
		

		tree0->Fill();	
		events_saved++;	
		float a1,a2,a3,a4,a5,a6,a7,a8,total;
		a1=0.15;
		a2=0.32;	
		a3=0.39;
		a4=0.095;	
		a5=0.024;	
		a6=0.011;	
		a7=0.0026;
		a8=0.00055;
		if  (region.CompareTo("mu")==0) total=180000;
		if  (region.CompareTo("el")==0) total=91000;

			
		if (file_tag.CompareTo("DYJetstoLL_madgraph")==0) if (events_saved>=total*a1) break;
		if (file_tag.CompareTo("DYJetstoLL_HT100to200")==0) if (events_saved>=total*a2) break;
		if (file_tag.CompareTo("DYJetstoLL_HT200to400")==0) if (events_saved>=total*a3) break;
		if (file_tag.CompareTo("DYJetstoLL_HT400to600")==0) if (events_saved>=total*a4) break;
		if (file_tag.CompareTo("DYJetstoLL_HT600to800")==0) if (events_saved>=total*a5) break;
		if (file_tag.CompareTo("DYJetstoLL_HT800to1200")==0) if (events_saved>=total*a6) break;
		if (file_tag.CompareTo("DYJetstoLL_HT1200to2500")==0) if (events_saved>=total*a7) break;
		if (file_tag.CompareTo("DYJetstoLL_HT2500toInf")==0) if (events_saved>=total*a8) break;

		if (file_tag_str.find("EWK")!=std::string::npos) if (region.CompareTo("mu")==0) if (events_saved>=total) break;
		if (file_tag_str.find("EWK")!=std::string::npos) if (region.CompareTo("el")==0) if (events_saved>=total) break;

	}  

	file.Write();
	file.Close();

}

