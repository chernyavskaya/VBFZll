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

fdy = ROOT.TFile.Open("/afs/cern.ch/work/n/nchernya/VBFZll/plotter/output_txt12/DYJetstoLL_amc_%s_QCDScalenom_JESnom_v25_bdt_alldata4_qglweightsnorm_cutbased2_reminiaod.root"%channel)
fewk = ROOT.TFile.Open("/afs/cern.ch/work/n/nchernya/VBFZll/plotter/output_txt12/EWK_LLJJ_%s_QCDScalenom_JESnom_v25_bdt_alldata4_qglweightsnorm_cutbased2_reminiaod.root"%channel)

hist_names = ['hMqq_bdt', 'hEtaQQ_bdt', 'hqgl_bdt', 'hqgl2_bdt','hqq_pt_bdt', 'hRptHard_bdt','hZll_zstar_bdt']
#hist_names = ['hMqq_bdt', 'hqgl2_bdt','hqq_pt_bdt', 'hRptHard_bdt','hZll_zstar_bdt']

hist_list_dy=[]
hist_list_ewk=[]

for i,name in enumerate(hist_names):
	print name
	hist_dy = fdy.Get(hist_names[i])
	hist_ewk = fewk.Get(hist_names[i])
	hist_dy.SetLineColor(2)
	hist_ewk.SetLineColor(4)
	hist_dy.SetFillColor(2)
	hist_ewk.SetFillColor(38)
	if (name=='hMqq_bdt' or name=='hEtaQQ_bdt' or name=='hqq_pt_bdt') : 
		hist_dy.Rebin(4)
		hist_ewk.Rebin(4)
	hist_dy.SetLineWidth(2)
	hist_ewk.SetLineWidth(2)
	hist_ewk.SetFillStyle(1001)
	hist_dy.SetFillStyle(3554)
	hist_list_dy.append(hist_dy)
	hist_list_ewk.append(hist_ewk)

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



for num in range(0,len(hist_names)):
#	list_amc[num].Scale(1./list_amc[num].Integral())
#	list_mdg[num].Scale(1./list_mdg[num].Integral())
	hist_list_dy[num].Scale(1./hist_list_dy[num].Integral())
	hist_list_ewk[num].Scale(1./hist_list_ewk[num].Integral())	
	c = ROOT.TCanvas("c","c",600,600)
	xmin = hist_list_dy[num].GetXaxis().GetXmin()
	xmax = hist_list_dy[num].GetXaxis().GetXmax()
	frame = ROOT.TH1F("frame","",1,xmin,xmax)
	frame.SetStats(0)
	frame.GetYaxis().SetTitle("Normalized Events")
	frame.GetYaxis().SetLabelSize(0.04)
	frame.GetXaxis().SetTitle(hist_list_dy[num].GetXaxis().GetTitle())
	frame.GetYaxis().SetRangeUser(0.,hist_list_ewk[num].GetMaximum()*1.3)
	frame.Draw()
	hist_list_ewk[num].Draw("HISTsame")
	hist_list_dy[num].Draw("HISTsame")

	leg = ROOT.TLegend(0.7,0.75,0.9,0.9)
	leg.AddEntry(hist_list_ewk[num],"Signal" ,"F")
	leg.AddEntry(hist_list_dy[num],"Background" ,"F")
	leg.SetFillStyle(-1)
	leg.SetBorderSize(0)
	leg.SetTextFont(42)
	leg.SetTextSize(0.03)
	leg.Draw()
	c.SaveAs("plots_after_bdt/plot_%s_%s.png" %(hist_names[num],channel) )





