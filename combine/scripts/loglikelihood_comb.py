import ROOT
from ROOT import gROOT
from ROOT import gStyle

gROOT.SetBatch(True)
gROOT.ProcessLineSync(".x /afs/cern.ch/work/n/nchernya/Hbb/setTDRStyle.C")
gROOT.ForceStyle()
gStyle.SetPadTopMargin(0.06)
gStyle.SetPadRightMargin(0.04)
gStyle.SetPadLeftMargin(0.15)


xmin=-0.01
xmax=0.01
#ymax=5
ymax=90
#end='FullStepSize'
end='_ForSign2'

path='/afs/cern.ch/work/n/nchernya/VBFZll/combine/run_combine_alldata5_qgl05/'
#f = ROOT.TFile.Open(path+"higgsCombinemuobs%s.MultiDimFit.mH125.root"%end)
f = ROOT.TFile.Open(path+"higgsCombinemu%sObs.MultiDimFit.mH125.root"%end)
tree = f.Get("limit")
gr = ROOT.TGraph()
for entry in range(1,tree.GetEntries()):
	tree.GetEntry(entry)
	gr.SetPoint(entry,tree.r,2*tree.deltaNLL)
	if tree.r==0.0005 : print tree.r, 2*tree.deltaNLL
#	print tree.r, 2*tree.deltaNLL

gr.SetLineColor(ROOT.kGreen+1)
gr.GetYaxis().SetTitle("-2 #Delta ln L")
gr.GetYaxis().SetRangeUser(0.,ymax)
gr.GetXaxis().SetLimits(xmin,xmax)
gr.GetXaxis().SetTitle("#mu")
gr.SetLineStyle(1)
gr.SetLineWidth(2)


#f2= ROOT.TFile.Open(path+"higgsCombinemuexp%s.MultiDimFit.mH125.root"%end)
f2= ROOT.TFile.Open(path+"higgsCombinemu%s.MultiDimFit.mH125.root"%end)
tree2 = f2.Get("limit")
gr2 = ROOT.TGraph()
for entry in range(1,tree2.GetEntries()):
	tree2.GetEntry(entry)
	gr2.SetPoint(entry,tree2.r,2*tree2.deltaNLL)
#	if entry==1: print tree2.r, 2*tree2.deltaNLL
	if tree2.r==0.0005 : print tree2.r, 2*tree2.deltaNLL
gr2.SetLineColor(ROOT.kGreen+1)
gr2.GetYaxis().SetTitle("-2 #Delta ln L")
ymax=7
gr2.GetYaxis().SetRangeUser(0.,ymax)
gr2.GetXaxis().SetLimits(xmin,xmax)
gr2.GetXaxis().SetTitle("#mu")
gr2.SetLineStyle(7)
gr2.SetLineWidth(2)

#f3= ROOT.TFile.Open(path+"higgsCombineelobs%s.MultiDimFit.mH125.root"%end)
f3= ROOT.TFile.Open(path+"higgsCombineel%sObs.MultiDimFit.mH125.root"%end)
tree3 = f3.Get("limit")
gr3 = ROOT.TGraph()
for entry in range(1,tree3.GetEntries()):
	tree3.GetEntry(entry)
	gr3.SetPoint(entry,tree3.r,2*tree3.deltaNLL)
#	if entry==1: print tree3.r, 2*tree3.deltaNLL
	if tree3.r==0.0005 : print tree3.r, 2*tree3.deltaNLL
gr3.SetLineColor(ROOT.kBlue+1)
gr3.SetLineStyle(1)
gr3.SetLineWidth(2)

#f4= ROOT.TFile.Open(path+"higgsCombineelexp%s.MultiDimFit.mH125.root"%end)
f4= ROOT.TFile.Open(path+"higgsCombineel%s.MultiDimFit.mH125.root"%end)
tree4 = f4.Get("limit")
gr4 = ROOT.TGraph()
for entry in range(1,tree4.GetEntries()):
	tree4.GetEntry(entry)
	gr4.SetPoint(entry,tree4.r,2*tree4.deltaNLL)
#	if entry==1: print tree4.r, 2*tree4.deltaNLL
	if tree4.r==0.0005 : print tree4.r, 2*tree4.deltaNLL
#	print tree4.r, 2*tree4.deltaNLL
gr4.SetLineColor(ROOT.kBlue+1)
gr4.SetLineStyle(7)
gr4.SetLineWidth(2)

#f5= ROOT.TFile.Open(path+"higgsCombineelmus%s.MultiDimFit.mH125.root"%end)
f5= ROOT.TFile.Open(path+"higgsCombineelmu%sObs.MultiDimFit.mH125.root"%end)
tree5 = f5.Get("limit")
gr5 = ROOT.TGraph()
for entry in range(1,tree5.GetEntries()):
	tree5.GetEntry(entry)
	gr5.SetPoint(entry,tree5.r,2*tree5.deltaNLL)
#	if entry==1: print tree5.r, 2*tree5.deltaNLL
	if tree5.r==0.0005 : print tree5.r, 2*tree5.deltaNLL
gr5.SetLineColor(ROOT.kRed)
gr5.SetLineStyle(1)
gr5.SetLineWidth(2)

#f6= ROOT.TFile.Open(path+"higgsCombineelmuexp%s.MultiDimFit.mH125.root"%end)
f6= ROOT.TFile.Open(path+"higgsCombineelmu%s.MultiDimFit.mH125.root"%end)
tree6 = f6.Get("limit")
gr6 = ROOT.TGraph()
for entry in range(1,tree6.GetEntries()):
	tree6.GetEntry(entry)
	gr6.SetPoint(entry,tree6.r,2*tree6.deltaNLL)
#	if entry==1: print tree6.r, 2*tree6.deltaNLL
	if tree6.r==0.0005 : print tree6.r, 2*tree6.deltaNLL
	print tree6.r, 2*tree6.deltaNLL
gr6.SetLineColor(ROOT.kRed)
gr6.SetLineStyle(7)
gr6.SetLineWidth(2)




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


gr.GetXaxis().SetLabelSize(0.04);
gr.GetYaxis().SetLabelSize(0.04);
gr.Draw("AC")
#gr2.Draw("AC")
gr2.Draw("Csame")
gr3.Draw("Csame")
gr4.Draw("Csame")
gr5.Draw("Csame")
gr6.Draw("Csame")
pCMS3 = ROOT.TPaveText(0.48,0.87,.65,.9,"NDC")
#pCMS3 = ROOT.TPaveText(0.73,0.87,.85,.9,"NDC")
pCMS3.SetTextFont(42)
pCMS3.SetTextSize(0.03)
pCMS3.SetFillStyle(1001)
#pCMS3.SetFillColor(ROOT.kWhite)
pCMS3.SetFillStyle(-1)
pCMS3.SetBorderSize(0)
pCMS3.AddText("EWK Zjj")
pCMS3.Draw()
leg = ROOT.TLegend(0.4,0.63,0.75,0.85)
#leg = ROOT.TLegend(0.65,0.63,0.85,0.85)
leg.SetFillStyle(-1)
leg.SetBorderSize(0)
leg.AddEntry(gr,"Dimuon Observed","L")
leg.AddEntry(gr2,"Dimuon Expected","L")
leg.AddEntry(gr3,"Dielectron Observed","L")
leg.AddEntry(gr4,"Dielectron Expected","L")
leg.AddEntry(gr5,"Combined Observed","L")
leg.AddEntry(gr6,"Combined Expected","L")

#observed
#leg.AddEntry(gr6,"8TeV","L")
#leg.AddEntry(gr,"13TeV","L")
#leg.AddEntry(gr3,"8TeV + 13TeV","L")

#expected
#leg.AddEntry(gr5,"8TeV","L")
#leg.AddEntry(gr2,"13TeV","L")
#leg.AddEntry(gr4,"8TeV + 13TeV","L")

leg.SetFillStyle(1001)
leg.SetFillColor(ROOT.kWhite)
leg.SetBorderSize(0)
leg.SetTextFont(42)
leg.SetTextSize(0.03)

line  = ROOT.TLine(xmin,1.,xmax,1.)
line.SetLineColor(1)
line.SetLineStyle(9)
line.SetLineWidth(1)
line.Draw()
pCMS4 = ROOT.TPaveText(0.203,0.76,.24,.78,"NDC")
pCMS4.SetTextFont(42)
pCMS4.SetTextSize(0.025)
pCMS4.SetFillStyle(-1)
pCMS4.SetBorderSize(0)
pCMS4.AddText("95% CL")
pCMS4.Draw()
line2  = ROOT.TLine(xmin,3.84159,xmax,3.84159)
line2.SetLineColor(1)
line2.SetLineStyle(9)
line2.SetLineWidth(1)
line2.Draw()
pCMS5 = ROOT.TPaveText(0.203,0.3,.24,.32,"NDC")
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

pCMS12 = ROOT.TPaveText(left+0.1,1.-top*1.13,0.57,1.,"NDC")
pCMS12.SetTextFont(52)
pCMS12.SetTextSize(top*0.73)
pCMS12.SetTextAlign(12)
pCMS12.SetFillStyle(-1)
pCMS12.SetBorderSize(0)
pCMS12.AddText("Preliminary")
pCMS12.Draw()


pCMS2 = ROOT.TPaveText(0.5,1.-top,1.-right*0.5,1.,"NDC")
pCMS2.SetTextFont(42)
pCMS2.SetTextSize(top*0.55)
pCMS2.SetTextAlign(32)
pCMS2.SetFillStyle(-1)
pCMS2.SetBorderSize(0)
pCMS2.AddText("L = 35.9 fb^{-1} (13 TeV)")
pCMS2.Draw()


leg.Draw()

c.SaveAs("plots/loglike_combined_%s_Obs.pdf"%end)
c.SaveAs("plots/loglike_combined_%s_Obs.png"%end)


