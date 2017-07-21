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
void efficiency_roc_Nm1(){
	gROOT->ProcessLine(".x /afs/cern.ch/work/n/nchernya/setTDRStyle.C");
	const int n_variables =8  ;
//	TString variables_names[n_variables]={"Mqq", "DeltaEtaQQ", "Jet2q_pt" ,"Jet1q_pt", "Jet1q_leadTrackPt", "Jet2q_leadTrackPt", "axis2_jet1", "axis2_jet2", "qq_pt", "RptHard", "Zll_pt", "Zll_zstar","all"}; 
//	TString variables_names[n_variables]={"Mqq", "DeltaEtaQQ", "Jet2q_pt" ,"axis2_jet1", "axis2_jet2", "qq_pt", "RptHard", "Zll_zstar","all"}; 
//	TString variables_names[n_variables]={"Mqq", "DeltaEtaQQ","axis2_jet1", "axis2_jet2", "qq_pt", "RptHard", "Zll_zstar","all"}; 
	TString variables_names[n_variables]={"Mqq", "DeltaEtaQQ","Jet2q_pt","axis2_jet1", "axis2_jet2", "qq_pt",  "Zll_zstar","all"}; 
//	TString variables_names[n_variables]={"Mqq", "DeltaEtaQQ", "axis2_jet1", "axis2_jet2",  "Zll_zstar","all"}; 
	TString file_names[n_variables];
	TString file_names_single[n_variables];
	Double_t frame2_axisx[n_variables];
	for (int i=0;i<n_variables;i++){
		frame2_axisx[i] = 0.5+i;
	}

	TString end = "axis2jet2q";
	Double_t hist_integrals[n_variables];
	Double_t hist_integrals_single[n_variables];
	for (int i=0;i<n_variables;i++){
		file_names[i].Prepend("/mnt/t3nfs01/data01/shome/nchernya/VBFZll/main_mva/output/v25/TMVA_main_v25_Nm1_"+variables_names[i]);
		file_names_single[i] = file_names[i];
		file_names[i].Append("elmu"+end+".root");
		file_names_single[i].Append("elmu"+end+".root");
	}

	for (int current_file=0;current_file<n_variables;current_file++){
//	for (int current_file=n_variables-1;current_file>=0;current_file--){
		TFile *file = new TFile(file_names[current_file]);
		file->cd("Method_BDT/BDTG");
		file->ls();
	//	TH1D *hist = (TH1D*) MVA_BDTG_rejBvsS->Clone(); 
		TH1D *hist = (TH1D*)file->Get("Method_BDT/BDTG/MVA_BDTG_rejBvsS"); 
		hist_integrals[current_file] = hist->Integral("width");
		cout<<hist_integrals[current_file]<<endl;
		TFile *file_single = new TFile(file_names_single[current_file]);
		file_single->cd("Method_BDT/BDTG");
		file_single->ls();
	//	TH1D *hist_single = (TH1D*) MVA_BDTG_rejBvsS->Clone(); 
//		TH1D *hist_single = (TH1D*)file_single->Get("MVA_BDTG_rejBvsS"); 
		TH1D *hist_single = (TH1D*)file_single->Get("Method_BDT/BDTG/MVA_BDTG_rejBvsS"); 
		hist_integrals_single[current_file] = hist_single->Integral("width");
	}


		TCanvas *c2 = new TCanvas();
		c2->SetBottomMargin(.12);
		c2->cd();
		TH1F *frame2 = new TH1F("frame2","",n_variables,0.,n_variables);
		frame2->SetMinimum(0.75);
      frame2->SetMaximum(.84);
      frame2->SetMinimum(hist_integrals[n_variables-1]*0.9);
      frame2->SetMaximum(hist_integrals[n_variables-1]*1.05);
      frame2->GetYaxis()->SetTitleOffset(1.4);
      frame2->GetXaxis()->SetTitleOffset(1.);
      frame2->SetStats(0);
      frame2->SetTitleFont(42,"x");
		frame2->SetTitleFont(42,"y");
      frame2->SetTitleSize(0.05, "XYZ");
		frame2->SetYTitle("ROC Area");
		frame2->SetXTitle("");	
		frame2->GetXaxis()->SetLabelSize(0.05);
		for (int i=0;i<n_variables;i++){
			frame2->GetXaxis()->SetBinLabel(i+1,variables_names[i]);
		}
		frame2->Draw();
		TGraph *gr = new TGraph(n_variables,frame2_axisx,hist_integrals);
		gr->SetMarkerStyle(20);
		gr->SetLineWidth(2);
		gr->Draw("PLsame");
		TLegend *leg = new TLegend(0.7,0.86,0.9,0.9);
		leg->SetFillColor(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.04);
	//	leg->AddEntry(gr,"Dimuon","PL");
		leg->AddEntry(gr,"Dilepton","PL");
		leg->Draw("same");
		TLine *line2 = new TLine(0.,hist_integrals[n_variables-1],(Double_t)n_variables,  hist_integrals[n_variables-1]);
		line2->SetLineStyle(2);
		line2->SetLineColor(2);
		line2->SetLineWidth(2);
		line2->Draw("Lsame");
		c2->Print("plots/v25/efficiency_v25_mu_"+end+".png");

		TCanvas *c3 = new TCanvas();
		c3->SetBottomMargin(.12);
		c3->cd();
      frame2->SetMinimum(hist_integrals_single[n_variables-1]*0.9);
//		frame2->SetMinimum(0.84);
  //    frame2->SetMaximum(.887);
      frame2->SetMaximum(hist_integrals_single[n_variables-1]*1.05);
		frame2->Draw();
		TGraph *gr_single = new TGraph(n_variables,frame2_axisx,hist_integrals_single);
		gr_single->SetMarkerStyle(25);
		gr_single->SetLineColor(kBlue);
		gr_single->SetMarkerColor(kBlue);
		gr_single->SetLineWidth(2);
		gr_single->Draw("PLsame");

		TLegend *leg1 = new TLegend(0.7,0.86,0.9,0.9);
		leg1->SetFillColor(0);
		leg1->SetBorderSize(0);
		leg1->SetTextSize(0.04);
		leg1->AddEntry(gr_single,"Dielectron","PL");
		leg1->Draw("same");
		TLine *line3 = new TLine(0.,hist_integrals_single[n_variables-1],(Double_t)n_variables, hist_integrals_single[n_variables-1]);
		line3->SetLineStyle(2);
		line3->SetLineColor(2);
		line3->SetLineWidth(2);
		line3->Draw("Lsame");
	
		c3->Print("plots/v25/efficiency_v25_el_"+end+".png");

}
