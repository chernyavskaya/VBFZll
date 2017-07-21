import ROOT
from ROOT import gROOT
from ROOT import gStyle
import sys

def getPt(item):
	return item[0]
def getEta(item):
	return item[1]

channel=sys.argv[1]

gROOT.SetBatch(True)
gROOT.ProcessLineSync(".x /afs/cern.ch/work/n/nchernya/Hbb/setTDRStyle.C")
gROOT.ForceStyle()
gStyle.SetPadTopMargin(0.06)
gStyle.SetPadRightMargin(0.04)
gStyle.SetPadLeftMargin(0.15)

#f = ROOT.TFile.Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat//store/user/nchernya/VBFZll/systematics/v25/interference_%s_v25_systematics_alldata4_qgl.root"%channel)
f = ROOT.TFile.Open("/afs/cern.ch/work/n/nchernya/VBFZll/combine/inputs/root/v25_alldata5/EWKzjj_v25_systematics_qgl_syst_interNew.root")
f2 = ROOT.TFile.Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat//store/user/nchernya/VBFZll/systematics/v25/LLJJ_EWK_5f_LO_13TeV_%s_v25_systematics_alldata4_qgl2.root"%channel)

#hist_names = ['atanhBDT_%s_interference'%channel,'atanhBDT_%s_interference_CMS_ewkzjj_int_shapeUp'%channel,'atanhBDT_%s_interference_CMS_ewkzjj_int_shapeDown'%channel]
hist_names = ['atanhBDT_%s_EWKinterference'%channel,'atanhBDT_%s_EWKinterference_CMS_ewkzjj_int_shapeUp'%channel,'atanhBDT_%s_EWKinterference_CMS_ewkzjj_int_shapeDown'%channel]
hist_names2 = ['atanhBDT_%s_LLJJ_EWK_5f_LO_13TeV'%channel]
list = []

for i,name in enumerate(hist_names):
	f.ls()
	hist = f.Get(hist_names[i])
	hist.SetLineWidth(2)
	hist.Rebin(30)
	list.append(hist)
list[0].SetLineColor(ROOT.kBlue)
list[1].SetLineColor(ROOT.kRed)
list[2].SetLineColor(ROOT.kGreen)
histInter = f2.Get(hist_names2[0])
histInter.SetLineWidth(2)
histInter.Rebin(30)
histInter.SetLineColor(ROOT.kViolet)


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
pCMS2.AddText("(13 TeV)")
	
pName = ROOT.TPaveText(left+0.01,1.-top*2,0.3,.9,"NDC")
pName.SetTextFont(42)
pName.SetTextSize(top*0.55)
pName.SetTextAlign(12)
pName.SetFillStyle(-1)
pName.SetBorderSize(0)
pName.AddText("%s channel"%channel)


if (1>0):
	c = ROOT.TCanvas("c","c",900,900)
	c.cd()
#	c.SetBottomMargin(0.3)
	xmin = list[0].GetXaxis().GetXmin()
	xmax = list[0].GetXaxis().GetXmax()
	frame = ROOT.TH1F("frame","",1,xmin,xmax)
	frame.SetStats(0)
	frame.GetXaxis().SetTitle("BDT'")
	frame.GetYaxis().SetTitle("Events")
	frame.GetYaxis().SetLabelSize(0.04)
	frame.GetYaxis().SetRangeUser(0.,list[1].GetMaximum()*1.15)
	frame.Draw()
	list[0].Draw("HISTsameE")
	list[1].Draw("HISTsameE")
	list[2].Draw("HISTsameE")
#	histInter.Draw("HISTsameE")
	pCMS1.Draw()
	pCMS12.Draw()
	pCMS2.Draw()

	leg = ROOT.TLegend(0.7,0.75,0.9,0.9)
	leg.AddEntry(list[0],"Nominal" ,"L")
	leg.AddEntry(list[1],"Up variation" ,"L")
	leg.AddEntry(list[2],"Down variation" ,"L")
	leg.SetFillStyle(-1)
	leg.SetBorderSize(0)
	leg.SetTextFont(42)
	leg.SetTextSize(0.03)
	leg.Draw()
	pName.Draw()

	ROOT.gPad.Update()
	ROOT.gPad.RedrawAxis()
	c.SaveAs("plots/new_interference_with_unc_shape_unc_%s_bdt.pdf"%channel)





