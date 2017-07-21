import ROOT
from ROOT import gROOT,TLine,TF1,TMath,TChain, TH1F
from ROOT import gStyle
import sys
import math


channel = sys.argv[1]

if channel=="mu" : channel_name = "Dimuon"
if channel=="el" : channel_name = "Dielectron"

gROOT.SetBatch(True)
gROOT.ProcessLineSync(".x /afs/cern.ch/work/n/nchernya/Hbb/setTDRStyle.C")
gROOT.ForceStyle()
gStyle.SetPadTopMargin(0.06)
gStyle.SetPadRightMargin(0.04)
gStyle.SetPadLeftMargin(0.15)

f = ROOT.TFile.Open("/afs/cern.ch/work/n/nchernya/VBFZll/plotter/output_txt14/EWK_LLJJ_%s_QCDScalenom_JESnom_v25_bdt_alldata4_qglweightsnorm_vetoeff_reminiaod.root"%channel)


path='/afs/cern.ch/user/n/nchernya/eos_mount/0/cms/store/group/phys_higgs/vbfHbb/V25_passall/EWK_LLJJ_MLL-50_MJJ-120_13TeV-madgraph-pythia8/VHBB_HEPPY_V25/170401_193326/0000/tree*.root'

	
hist_sel = f.Get("hPVs")
hist_sel.Scale(1./hist_sel.Integral())

chain = TChain("tree")
chain.Add(path)
hnPVs = TH1F("hPVs_all","",50,0,50)
hnPVs.GetXaxis().SetTitle("nPVs")
chain.Draw("nPVs>>hPVs_all","puWeight*genWeight/TMath::Abs(genWeight)")
hnPVs.Scale(1./hnPVs.Integral())

hnPVs.SetLineColor(2)
hnPVs.SetLineWidth(2)
hist_sel.SetLineColor(ROOT.kBlue)
hist_sel.SetLineWidth(2)
hist_sel.SetLineStyle(7)


hratio = hist_sel.Clone("new")
hratio.Divide(hnPVs)



right = gStyle.GetPadRightMargin()
top = gStyle.GetPadTopMargin()
left = gStyle.GetPadLeftMargin()
bottom =  gStyle.GetPadBottomMargin()

pCMS1 = ROOT.TPaveText(left,1.-top,0.4,1.,"NDC")
pCMS1.SetTextFont(62)
pCMS1.SetTextSize(top*0.75)
pCMS1.SetTextAlign(12)
pCMS1.SetFillStyle(-1)
pCMS1.SetBorderSize(0)
pCMS1.AddText("CMS")

pCMS12 = ROOT.TPaveText(left+0.1,1.-top*1.13,0.57,1.,"NDC")
pCMS12.SetTextFont(52)
pCMS12.SetTextSize(top*0.73)
pCMS12.SetTextAlign(12)
pCMS12.SetFillStyle(-1)
pCMS12.SetBorderSize(0)
pCMS12.AddText("Work in progress")

pCMS11 = ROOT.TPaveText(left,1.-top*2,0.4,1.-top*1.5,"NDC")
pCMS11.SetTextFont(42)
pCMS11.SetTextSize(top*0.75)
pCMS11.SetTextAlign(12)
pCMS11.SetFillStyle(-1)
pCMS11.SetBorderSize(0)
pCMS11.AddText(channel_name)


pCMS2 = ROOT.TPaveText(0.5,1.-top,1.-right*0.5,1.,"NDC")
pCMS2.SetTextFont(42)
pCMS2.SetTextSize(top*0.75)
pCMS2.SetTextAlign(32)
pCMS2.SetFillStyle(-1)
pCMS2.SetBorderSize(0)
pCMS2.AddText("35.9 fb^{-1} (13 TeV)")

if (1>0):	
	c = ROOT.TCanvas("c","c",900,900)
	c.cd()
	c.SetBottomMargin(0.3)
	xmin = hist_sel.GetXaxis().GetXmin()*0.9
	xmax = hist_sel.GetXaxis().GetXmax()*1.1
	frame = ROOT.TH1F("frame","",1,xmin,xmax)
	frame.SetStats(0)
	frame.GetYaxis().SetTitle("Normalized Events")
	frame.GetYaxis().SetLabelSize(0.04)
	frame.GetXaxis().SetLabelSize(0)
	frame.GetXaxis().SetTitleOffset(0.91);
#	frame.GetXaxis().SetTitle(list_hist[num].GetXaxis().GetTitle())
	y_min= hist_sel.GetMinimum()*0.9
	y_max = hist_sel.GetMaximum()*1.2
	frame.GetYaxis().SetRangeUser(y_min,y_max)
	frame.Draw()
	hnPVs.Draw("HISTsame")
	hist_sel.Draw("HIStsame")
	pCMS1.Draw()
	pCMS11.Draw()
	pCMS12.Draw()
	pCMS2.Draw()
	leg = ROOT.TLegend(0.6,0.75,0.9,0.9)
	leg.AddEntry(hist_sel,"after preselection" ,"L")
	leg.AddEntry(hnPVs,"no preselection" ,"L")
	leg.SetFillStyle(-1)
	leg.SetBorderSize(0)
	leg.SetTextFont(42)
	leg.SetTextSize(0.03)
	leg.Draw()


	pad2 = ROOT.TPad("pad2", "pad2", 0., 0., 1., 1.)
	pad2.SetTopMargin(0.73)
#	pad2.SetBottomMargin(0.)
#	pad2.SetRightMargin(0.)
	pad2.SetFillColor(0)
	pad2.SetFillStyle(0)
	pad2.Draw()
	pad2.cd()
	frame2 = ROOT.TH1F("frame2","",1,xmin,xmax)
	frame2.SetMinimum(0.)	
	frame2.SetMaximum(2) 
	frame2.GetYaxis().SetLabelSize(0.02)
	frame2.GetXaxis().SetLabelSize(0.04)
	frame2.GetYaxis().SetTitleSize(0.04)
	frame2.GetXaxis().SetTitle(hist_sel.GetXaxis().GetTitle())
	frame2.SetStats(0)
	frame2.GetYaxis().SetTitle("Presel/All")	
	frame2.Draw()

	hratio.Draw("PEsame")
		
	line2 = TLine(xmin,1,xmax,1)
	line2.SetLineStyle(2)
	line2.SetLineColor(1)
	line2.SetLineWidth(2)
	line2.Draw("same")
	
	
	ROOT.gPad.Update()
	ROOT.gPad.RedrawAxis()
	c.SaveAs("pu_dep/plot_pudependence_%s.png" %(channel ))

