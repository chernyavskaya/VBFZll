#include <TH1D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TLorentzVector.h"
#include <iomanip>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "TFile.h"

void efficiency_groups(){
	const int n_groups = 6;
	TString file_names[n_groups] = {"1","2","3","4","5","6"};
	TString legend_names[n_groups] = {"qq system", "soft activity", "btag", "angular dynamics","bb system", "qb system"};
	Double_t frame2_axisx[n_groups] = {.5,1.5,2.5,3.5,4.5,5.5};
	Double_t hist_integrals[n_groups];
	for (int i=0;i<n_groups;i++){
		file_names[i].Prepend("../output/TMVA_main_");
		//file_names[i].Append("_step_qcd_300to500_double.root");
		file_names[i].Append("_step_qcd_300to500_single.root");
	}
TLegend *leg = new TLegend(0.12,0.17,0.45,0.4);
leg->SetBorderSize(0);
leg->SetTextSize(0.04);


		TLatex* tex = new TLatex(0.90,0.92,"13 TeV, PU = 20, bx = 25 ns");
      tex->SetNDC();
		tex->SetTextAlign(35);
      tex->SetTextFont(42);
      tex->SetTextSize(0.04);
      tex->SetLineWidth(2);
      TLatex *tex1 = new TLatex(0.13,0.75,"CMS");
      tex1->SetNDC();
      tex1->SetTextAlign(20);
      tex1->SetTextFont(61);
      tex1->SetTextSize(0.06);
      tex1->SetLineWidth(2);
      TLatex* tex2 = new TLatex(0.22,0.68,"Work in progress");
      tex2->SetNDC();
      tex2->SetTextAlign(20);
      tex2->SetTextFont(52);
      tex2->SetTextSize(0.04);
      tex2->SetLineWidth(2);
	//	TLatex* tex_file = new TLatex(0.35,0.92,"Spring15, DoubleBtag");
		TLatex* tex_file = new TLatex(0.35,0.92,"Spring15, SingleBtag");
      tex_file->SetNDC();
		tex_file->SetTextAlign(35);
      tex_file->SetTextFont(42);
      tex_file->SetTextSize(0.04);
      tex_file->SetLineWidth(2);	
		TCanvas *c1 = new TCanvas();
		c1->SetBottomMargin(.12);
		c1->cd();
		TH1F *frame = new TH1F("frame","",1,0.,1.);
		frame->SetMinimum(0.);
      frame->SetMaximum(1.05);
      frame->GetYaxis()->SetTitleOffset(0.9);
      frame->GetXaxis()->SetTitleOffset(0.91);
      frame->SetStats(0);
      frame->SetTitleFont(42,"x");
		frame->SetTitleFont(42,"y");
      frame->SetTitleSize(0.05, "XYZ");
		frame->SetYTitle("Background rejection");
		frame->SetXTitle("Signal efficiency");
		frame->Draw();
	for (int current_file=n_groups-1;current_file>=0;current_file--){
		TFile *file = new TFile(file_names[current_file]);
		file->cd("Method_BDT/BDTG");
		file->ls();
		TH1D *hist = (TH1D*) MVA_BDTG_rejBvsS->Clone(); 
		hist_integrals[current_file] = hist->Integral("width");
		hist->SetDirectory(0);
		hist->SetTitle("");
		hist->SetLineWidth(3);
		if (current_file==4) {
			hist->SetLineColor(kGreen);
			hist->SetLineStyle(2);
		}
		if (current_file==0) {
			hist->SetLineColor(kCyan);
			hist->SetLineStyle(5);
		}
		if (current_file==1) {
			hist->SetLineColor(kBlue);
			hist->SetLineStyle(7);
		}
		if (current_file==2) {
			hist->SetLineColor(1);
			hist->SetLineStyle(9);
		}
		if (current_file==3) {
			hist->SetLineColor(kOrange+1);
			hist->SetLineStyle(10);
		}
		if (current_file==5) {
			hist->SetLineColor(kRed);
			hist->SetLineStyle(1);
		}
		hist->Draw("Lsame");
		leg->AddEntry(hist, legend_names[current_file]);
		tex->Draw();
		tex1->Draw();
		tex2->Draw();
		tex_file->Draw();
		//file->Close();
	}
	leg->Draw("same");
	c1->Print("plots/efficiency_groups_roc_single.png");


}
