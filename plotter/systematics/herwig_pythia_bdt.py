import ROOT
from ROOT import gROOT
from ROOT import gStyle
import sys

channel=sys.argv[1]
bdt=sys.argv[2]

if channel=="mu" : channel2="Dimuon"
if channel=="el" : channel2="Dielectron"



gROOT.SetBatch(True)
gROOT.ProcessLineSync(".x /afs/cern.ch/work/n/nchernya/Hbb/setTDRStyle.C")
gROOT.ForceStyle()
gStyle.SetPadTopMargin(0.06)
gStyle.SetPadRightMargin(0.04)
gStyle.SetPadLeftMargin(0.15)


right,top   = gStyle.GetPadRightMargin(),gStyle.GetPadTopMargin()
left,bottom = gStyle.GetPadLeftMargin(),gStyle.GetPadBottomMargin()
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
pCMS12.AddText("Simulation")



pCMS2 = ROOT.TPaveText(0.5,1.-top,1.-right*0.5,1.,"NDC")
pCMS2.SetTextFont(42)
pCMS2.SetTextSize(top*0.75)
pCMS2.SetTextAlign(32)
pCMS2.SetFillStyle(-1)
pCMS2.SetBorderSize(0)
pCMS2.AddText("35.9 fb^{-1} (13 TeV)")




file_names = ['EWK_LLJJ','EWK_LLJJ_herwig']

hist_list=[]

for i in range(0,len(file_names)):
	f = ROOT.TFile.Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat//store/user/nchernya/VBFZll/plotterOutput/v25/%s_%s_QCDScalenom_JESnom_v25_bdt_alldata4_qglweightsnorm_cutbased_reminiaod.root"%(file_names[i],channel))
	hist = f.Get(bdt)
	hist.Scale(35900.)
	if (bdt!="hbdt") : hist.Rebin(10)
	if (bdt=="hbdt") : hist.Rebin(2)
	hist_list.append(hist)


c = ROOT.TCanvas("c","c",900,900)
c.cd()
c.SetBottomMargin(0.15)
xmin = hist_list[0].GetXaxis().GetXmin()
xmax = hist_list[0].GetXaxis().GetXmax()
if (bdt!="hbdt") : xmax=3
frame = ROOT.TH1F("frame","",1,xmin,xmax)
frame.SetStats(0)
frame.GetXaxis().SetTitleOffset(0.91);
frame.GetYaxis().SetTitle("Events")
if (bdt=="hbdt") : frame.GetXaxis().SetTitle("BDT")
else : frame.GetXaxis().SetTitle("BDT'")
frame.GetYaxis().SetLabelSize(0.04)
frame.GetYaxis().SetRangeUser(0.,hist_list[0].GetMaximum()*1.4)
frame.Draw()
pCMS1.Draw()
pCMS12.Draw()
pCMS2.Draw()
pName = ROOT.TPaveText(left+0.01,1.-top*3,0.3,.9,"NDC")
pName.SetTextFont(42)
pName.SetTextSize(top*0.55)
pName.SetTextAlign(12)
pName.SetFillStyle(-1)
pName.SetBorderSize(0)
pName.AddText("%s"%channel2)
pName.Draw()

colors=[ROOT.kRed,ROOT.kBlack]
for i in range(0,len(file_names)):
	hist_list[i].SetLineColor(colors[i])
	hist_list[i].SetMarkerStyle(1)
	hist_list[i].SetLineWidth(2)
	hist_list[i].Draw("HSAME")
	
leg = ROOT.TLegend(0.7,0.75,0.9,0.9)
leg.AddEntry(hist_list[0],"Pythia" ,"L")
leg.AddEntry(hist_list[1],"Herwig" ,"L")
leg.SetFillStyle(-1)
leg.SetBorderSize(0)
leg.SetTextFont(42)
leg.SetTextSize(0.03)
leg.Draw()

c.SaveAs("plots_herwig/herwig_pythia_%s_%s.png"%(bdt,channel))


