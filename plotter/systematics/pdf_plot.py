import ROOT
from ROOT import gROOT
from ROOT import gStyle
import sys

channel=sys.argv[1]


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




file_names = ['EWK_LLJJ','TT','DYJetstoLL_amc_0J','DYJetstoLL_amc_1J','DYJetstoLL_amc_2J']

for i in range(0,len(file_names)):
	f = ROOT.TFile.Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat//store/user/nchernya/VBFZll/plotterOutput/v25/%s_%s_v25_pdfAcceptance_pdf_unc_reminiaod.root"%(file_names[i],channel))
	if (channel=='mu') and (file_names[i]=='DYJetstoLL_amc_2J'):
		f = ROOT.TFile.Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat//store/user/nchernya/VBFZll/plotterOutput/v25/%s_%s_v25_pdfAcceptance_pdf_unc2try_reminiaod.root"%(file_names[i],channel))
	hist = f.Get('acceptance')
	c = ROOT.TCanvas("c","c",900,900)
	c.cd()
	c.SetBottomMargin(0.15)
	hist.Draw("SAME")
	pCMS1.Draw()
	pCMS12.Draw()
	pCMS2.Draw()
	

	pMean = ROOT.TPaveText(0.7,0.8,0.9,.9,"NDC")
	pMean.SetTextFont(42)
	pMean.SetTextSize(top*0.55)
	pMean.SetTextAlign(32)
	pMean.SetFillStyle(-1)
	pMean.SetBorderSize(0)
	pRMS = ROOT.TPaveText(0.7,0.7,0.9,.8,"NDC")
	pRMS.SetTextFont(42)
	pRMS.SetTextSize(top*0.55)
	pRMS.SetTextAlign(32)
	pRMS.SetFillStyle(-1)
	pRMS.SetBorderSize(0)
	pRatio = ROOT.TPaveText(0.7,0.6,0.9,.7,"NDC")
	pRatio.SetTextFont(42)
	pRatio.SetTextSize(top*0.55)
	pRatio.SetTextAlign(32)
	pRatio.SetFillStyle(-1)
	pRatio.SetBorderSize(0)
	ratio = hist.GetRMS()/hist.GetMean()
	pMean.AddText('Mean = %0.3f'%hist.GetMean())
	pRMS.AddText('RMS = %0.3f'%hist.GetRMS())
	pRatio.AddText('RMS/Mean = %0.3f'%ratio)
	pMean.Draw()
	pRMS.Draw()
	pRatio.Draw()
	pName = ROOT.TPaveText(left+0.01,1.-top*3,0.3,.9,"NDC")
	pName.SetTextFont(42)
	pName.SetTextSize(top*0.55)
	pName.SetTextAlign(12)
	pName.SetFillStyle(-1)
	pName.SetBorderSize(0)
	pName.AddText("%s"%file_names[i])
	pName.AddText("%s channel"%channel)
	pName.Draw()
	c.SaveAs("plots_pdf/pdf_%s_%s.png"%(file_names[i],channel))
	c.SaveAs("plots_pdf/pdf_%s_%s.pdf"%(file_names[i],channel))


