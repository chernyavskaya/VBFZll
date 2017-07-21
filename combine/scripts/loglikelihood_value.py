import ROOT
from ROOT import gROOT
from ROOT import gStyle

gROOT.SetBatch(True)
gROOT.ProcessLineSync(".x /afs/cern.ch/work/n/nchernya/Hbb/setTDRStyle.C")
gROOT.ForceStyle()
gStyle.SetPadTopMargin(0.06)
gStyle.SetPadRightMargin(0.04)
gStyle.SetPadLeftMargin(0.15)

f = ROOT.TFile.Open("higgsCombinenewoption_nom200_2_15_7.MultiDimFit.mH125.root")
tree = f.Get("limit")
gr = ROOT.TGraph()
for entry in range(1,tree.GetEntries()):
	tree.GetEntry(entry)
	gr.SetPoint(entry,tree.r,2*tree.deltaNLL)
	print tree.r, 2*tree.deltaNLL

gr.SetLineColor(1)
gr.GetYaxis().SetTitle("-2 #Delta ln L")
ymax=7
gr.GetYaxis().SetRangeUser(0.,ymax)
gr.GetXaxis().SetLimits(-12.,9.)
gr.GetXaxis().SetTitle("#mu")
gr.SetLineStyle(1)
gr.SetLineWidth(2)


f2= ROOT.TFile.Open("higgsCombinenewoption_nom_exp.MultiDimFit.mH125.root")
tree2 = f2.Get("limit")
gr2 = ROOT.TGraph()
for entry in range(1,tree2.GetEntries()):
	tree2.GetEntry(entry)
	gr2.SetPoint(entry,tree2.r,2*tree2.deltaNLL)
	print tree2.r, 2*tree2.deltaNLL
gr2.SetLineColor(1)
gr2.SetLineStyle(7)
gr2.SetLineWidth(2)

c = ROOT.TCanvas("c","c",900,900)
#c.SetBottomMargin(.3);
##c.SetLeftMargin(.17);
#frame = ROOT.TH1F("frame","",60,-12.,9.);
#frame.SetStats(0);
#frame.GetXaxis().SetNdivisions(505);
#frame.GetXaxis().SetLabelSize(0.0);
#frame.SetYTitle("-2 #Delta ln L");
#frame.SetXTitle("#mu");
#frame.Draw()


gr.Draw("AC")
gr2.Draw("Csame")
pCMS3 = ROOT.TPaveText(0.5,0.8,.6,.9,"NDC")
pCMS3.SetTextFont(42)
pCMS3.SetTextSize(0.03)
pCMS3.SetFillStyle(-1)
pCMS3.SetBorderSize(0)
pCMS3.AddText("VBF H #rightarrow b#bar{b}")
pCMS3.Draw()
leg = ROOT.TLegend(0.47,0.73,0.63,0.82)
leg.AddEntry(gr,"Observed","L")
leg.AddEntry(gr2,"Expected","L")
leg.SetFillStyle(-1)
leg.SetBorderSize(0)
leg.SetTextFont(42)
leg.SetTextSize(0.03)
leg.Draw()

line  = ROOT.TLine(-12.,1.,9.,1.)
line.SetLineColor(1)
line.SetLineStyle(9)
line.SetLineWidth(1)
line.Draw()
pCMS4 = ROOT.TPaveText(0.203,0.59,.24,.61,"NDC")
pCMS4.SetTextFont(42)
pCMS4.SetTextSize(0.025)
pCMS4.SetFillStyle(-1)
pCMS4.SetBorderSize(0)
pCMS4.AddText("95% CL")
pCMS4.Draw()
line2  = ROOT.TLine(-12,3.84159,9,3.84159)
line2.SetLineColor(1)
line2.SetLineStyle(9)
line2.SetLineWidth(1)
line2.Draw()
pCMS5 = ROOT.TPaveText(0.203,0.26,.24,.28,"NDC")
pCMS5.SetTextFont(42)
pCMS5.SetTextSize(0.025)
pCMS5.SetFillStyle(-1)
pCMS5.SetBorderSize(0)
pCMS5.AddText("68% CL")
pCMS5.Draw()



right,top   = gStyle.GetPadRightMargin(),gStyle.GetPadTopMargin()
left,bottom = gStyle.GetPadLeftMargin(),gStyle.GetPadBottomMargin()
pCMS1 = ROOT.TPaveText(left,1.-top,0.4,1.,"NDC")
pCMS1.SetTextFont(62)
pCMS1.SetTextSize(top*0.75)
pCMS1.SetTextAlign(12)
pCMS1.SetFillStyle(-1)
pCMS1.SetBorderSize(0)
pCMS1.AddText("CMS")
pCMS1.Draw()
pCMS2 = ROOT.TPaveText(0.5,1.-top,1.-right*0.5,1.,"NDC")
pCMS2.SetTextFont(42)
pCMS2.SetTextSize(top*0.75)
pCMS2.SetTextAlign(32)
pCMS2.SetFillStyle(-1)
pCMS2.SetBorderSize(0)
pCMS2.AddText("L = 2.32 fb^{-1} (13 TeV)")
pCMS2.Draw()



c.SaveAs("loglike.pdf")
c.SaveAs("loglike.png")


