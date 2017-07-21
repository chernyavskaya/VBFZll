#include <TH1D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TLorentzVector.h"
#include <iomanip>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "TFile.h"

void efficiency_roc_groups(){
	const int n_groups = 6;
	TString file_names[n_groups] = {"1","2","3","4","5","6"};
	TString file_names_single[n_groups];
	TString legend_names[n_groups] = {"qq system", "soft activity", "btag", "angular dynamics","bb system", "qb system"};
	Double_t frame2_axisx[n_groups] = {.5,1.5,2.5,3.5,4.5,5.5};
	Double_t hist_integrals[n_groups];
	Double_t hist_integrals_single[n_groups];
	for (int i=0;i<n_groups;i++){
		file_names[i].Prepend("../output/TMVA_main_");
		file_names_single[i] = file_names[i];
		file_names[i].Append("_step_qcd_300to500_double.root");
		file_names_single[i].Append("_step_qcd_300to500_single.root");
	}
TLegend *leg = new TLegend(0.12,0.17,0.45,0.4);
leg->SetBorderSize(0);
leg->SetTextSize(0.04);

	for (int current_file=n_groups-1;current_file>=0;current_file--){
		TFile *file = new TFile(file_names[current_file]);
		file->cd("Method_BDT/BDTG");
		file->ls();
		TH1D *hist = (TH1D*) MVA_BDTG_rejBvsS->Clone(); 
		hist_integrals[current_file] = hist->Integral("width");
		TFile *file_single = new TFile(file_names_single[current_file]);
		file_single->cd("Method_BDT/BDTG");
		file_single->ls();
		TH1D *hist_single = (TH1D*) MVA_BDTG_rejBvsS->Clone(); 
		hist_integrals_single[current_file] = hist_single->Integral("width");
	}


		TCanvas *c2 = new TCanvas();
		c2->SetBottomMargin(.12);
		c2->cd();
		TH1F *frame2 = new TH1F("frame2","",n_groups,0.,6.);
		frame2->SetMinimum(0.6);
      frame2->SetMaximum(1.0);
      frame2->GetYaxis()->SetTitleOffset(0.9);
      frame2->GetXaxis()->SetTitleOffset(0.91);
      frame2->SetStats(0);
      frame2->SetTitleFont(42,"x");
		frame2->SetTitleFont(42,"y");
      frame2->SetTitleSize(0.05, "XYZ");
		frame2->SetYTitle("ROC Area");
		frame2->SetXTitle("");	
		frame2->GetXaxis()->SetLabelSize(0.05);
		for (int i=0;i<n_groups;i++){
			frame2->GetXaxis()->SetBinLabel(i+1,legend_names[i]);
		}
		frame2->Draw();
		TGraph *gr = new TGraph(n_groups,frame2_axisx,hist_integrals);
		gr->SetMarkerStyle(20);
		gr->SetLineWidth(2);
		gr->Draw("PLsame");
		TGraph *gr_single = new TGraph(n_groups,frame2_axisx,hist_integrals_single);
		gr_single->SetMarkerStyle(25);
		gr_single->SetLineColor(kBlue);
		gr_single->SetMarkerColor(kBlue);
		gr_single->SetLineWidth(2);
		gr_single->Draw("PLsame");

		TLegend *leg = new TLegend(0.6,0.15,0.9,0.4);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.04);
		leg->AddEntry(gr,"DoubleBtag","PL");
		leg->AddEntry(gr_single,"SingleBtag","PL");
		leg->Draw("same");
	
		c2->Print("plots/efficiency_groups_Nm1.png");

}
