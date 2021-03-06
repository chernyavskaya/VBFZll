import ROOT
from ROOT import TH1F, TFile
from ROOT import gROOT
from ROOT import gStyle
import sys

def getPt(item):
	return item[0]
def getEta(item):
	return item[1]

dy = sys.argv[1]

dy_choice=''
dy_model=''
if dy=='amc' : 
	dy_model = 'AMC@NLO'
	dy_choice = "amc_"
if dy=='mdg' : 
	dy_model = 'Madgraph'
	dy_choice = "HT_"


gROOT.SetBatch(True)
gROOT.ProcessLineSync(".x /afs/cern.ch/work/n/nchernya/Hbb/setTDRStyle.C")
gROOT.ForceStyle()
gStyle.SetPadTopMargin(0.06)
gStyle.SetPadRightMargin(0.04)
gStyle.SetPadLeftMargin(0.15)

path="dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat//store/user/nchernya/VBFZll/plotterOutput/v25/"
end="_QCDScalenom_JESnom_v25_bdt_alldata.root"

list_nom_data = []
list_denom_data = []
list_nom_dy = []
list_denom_dy = []
list_nom_ewk = []
list_denom_ewk = []
list_nom_ewk_herwig = []
list_denom_ewk_herwig = []

list_new_data = []
list_new_dy = []
list_new_ewk = []
list_new_ewk_herwig = []
list_new_dy_ewk = []
list_new_dy = []
list_new_dy_ewk_herwig = []

channel=["mu","el"]
hist_names_nom = ['hveto_jet3pt_nom', 'hveto_ht_nom', 'hveto_softht_nom', 'hveto_softpt_nom']
hist_names_denom = [ 'hveto_jet3pt_denom', 'hveto_ht_denom',  'hveto_softht_denom', 'hveto_softpt_denom']
axisx_names = ["Third jet p_{T} (GeV)","H_{T} (p_{T} > 20) (GeV)","Soft H_{T} (GeV)","Leading soft jet p_{T} (GeV)"]


for channel_name in channel:
	fdata = ROOT.TFile.Open(path+"SingleMuon_"+channel_name+end)
	if channel_name=="el" : fdata=ROOT.TFile.Open(path+"SingleElectron_"+channel_name+end)
	fdy = ROOT.TFile.Open(path+"DYJetstoLL_"+dy_choice+channel_name+end)
	fewk = ROOT.TFile.Open(path+"EWK_LLJJ_"+channel_name+end)
	fewk_herwig = ROOT.TFile.Open(path+"EWK_LLJJ_herwig_"+channel_name+end)
	for hist_name in hist_names_nom:
		hist = fdata.Get(hist_name)
		list_nom_data.append(hist)
		hist = fdy.Get(hist_name)
		list_nom_dy.append(hist)
		hist = fewk.Get(hist_name)
		list_nom_ewk.append(hist)
		hist = fewk_herwig.Get(hist_name)
		list_nom_ewk_herwig.append(hist)
	for hist_name in hist_names_denom:
		hist = fdata.Get(hist_name)
		list_denom_data.append(hist)
		hist = fdy.Get(hist_name)
		list_denom_dy.append(hist)
		hist = fewk.Get(hist_name)
		list_denom_ewk.append(hist)
		hist = fewk_herwig.Get(hist_name)
		list_denom_ewk_herwig.append(hist)



for i in range(0,len(hist_names_denom)):
	hist_new = list_denom_data[i].Clone("new_data") 
	hist_new.Add(list_nom_data[i],-1)
	hist_new.Divide(list_denom_data[i])
	list_new_data.append(hist_new)

	hist_new = list_denom_dy[i].Clone("new_dy")
	hist_new.Add(list_nom_dy[i],-1)
	hist_new.Divide(list_denom_dy[i])
	list_new_dy.append(hist_new)

	hist_new_nom = list_denom_dy[i].Clone("new_dy_nom")
	hist_new_nom.Add(list_nom_dy[i],-1)
	hist_new_nom.Add(list_denom_ewk[i])
	hist_new_nom.Add(list_nom_ewk[i],-1)
	hist_new_denom = list_denom_dy[i].Clone("new_dy_denom")
	hist_new_denom.Add(list_denom_ewk[i])
	hist_new_nom.Divide(hist_new_denom)
	list_new_dy_ewk.append(hist_new_nom)

	hist_new_nom = list_denom_dy[i].Clone("new_dy_nom")
	hist_new_nom.Add(list_nom_dy[i],-1)
	hist_new_nom.Add(list_denom_ewk_herwig[i])
	hist_new_nom.Add(list_nom_ewk_herwig[i],-1)
	hist_new_denom = list_denom_dy[i].Clone("new_dy_denom")
	hist_new_denom.Add(list_denom_ewk_herwig[i])
	hist_new_nom.Divide(hist_new_denom)
	list_new_dy_ewk_herwig.append(hist_new_nom)



right,top   = gStyle.GetPadRightMargin(),gStyle.GetPadTopMargin()
left,bottom = gStyle.GetPadLeftMargin(),gStyle.GetPadBottomMargin()
pCMS1 = ROOT.TPaveText(left,1.-top,0.4,1.,"NDC")
pCMS1.SetTextFont(62)
pCMS1.SetTextSize(top*0.75)
pCMS1.SetTextAlign(12)
pCMS1.SetFillStyle(-1)
pCMS1.SetBorderSize(0)
pCMS1.AddText("CMS")

pCMS12 = ROOT.TPaveText(left,1.-top*3.5,0.57,1.,"NDC")
pCMS12.SetTextFont(52)
pCMS12.SetTextSize(0.02)
pCMS12.SetTextAlign(12)
pCMS12.SetFillStyle(-1)
pCMS12.SetBorderSize(0)
pCMS12.AddText("ee + #mu#mu events, BDT > 0.92")
pCMS13 = ROOT.TPaveText(left,1.-top*4.6,0.57,1.,"NDC")
pCMS13.SetTextFont(52)
pCMS13.SetTextSize(0.02)
pCMS13.SetTextAlign(12)
pCMS13.SetFillStyle(-1)
pCMS13.SetBorderSize(0)
pCMS13.AddText("DY (%s)"%dy_model)



pCMS2 = ROOT.TPaveText(0.5,1.-top,1.-right*0.5,1.,"NDC")
pCMS2.SetTextFont(42)
pCMS2.SetTextSize(top*0.75)
pCMS2.SetTextAlign(32)
pCMS2.SetFillStyle(-1)
pCMS2.SetBorderSize(0)
pCMS2.AddText("35.9 fb^{-1} (13 TeV)")


for i,name in enumerate(hist_names_denom):
	c = ROOT.TCanvas("c","c",900,900)
	c.cd()
	c.SetBottomMargin(0.15)
	xmin = list_new_data[i].GetXaxis().GetXmin() -  list_new_data[i].GetXaxis().GetBinWidth(1)
	xmax = list_new_data[i].GetXaxis().GetXmax() +  list_new_data[i].GetXaxis().GetBinWidth(1)
	frame = ROOT.TH1F("frame","",1,xmin,xmax)
	frame.SetStats(0)
	frame.GetXaxis().SetTitleOffset(1.03)
	frame.GetXaxis().SetTitle(axisx_names[i])
	frame.GetYaxis().SetTitle("Gap fraction")
	frame.GetYaxis().SetLabelSize(0.04)
	frame.GetYaxis().SetRangeUser(list_new_data[i].GetMinimum()*0.7,1.1)
#	frame.GetYaxis().SetRangeUser(0.6,1.1)
	frame.Draw()
	list_new_dy_ewk[i].SetLineColor(ROOT.kViolet-5)
	list_new_dy_ewk[i].SetMarkerStyle(1)
	list_new_dy_ewk[i].SetLineWidth(2)
	hist_dy_ewk=list_new_dy_ewk[i].Clone("new1")
	list_new_dy_ewk[i].SetFillColor(ROOT.kViolet-9)
#	list_new_dy_ewk[i].SetFillStyle(4050)
	list_new_dy[i].SetLineColor(ROOT.kAzure+2)
	list_new_dy[i].SetLineWidth(2)
	list_new_dy[i].SetLineStyle(7)
	list_new_dy[i].SetMarkerStyle(1)
	hist_dy=list_new_dy[i].Clone("new2")
	list_new_dy[i].SetFillColor(ROOT.kAzure-9)
#	list_new_dy[i].SetFillStyle(3001)
	list_new_dy_ewk_herwig[i].SetLineColor(ROOT.kViolet+7)
	list_new_dy_ewk_herwig[i].SetLineWidth(2)
	list_new_dy_ewk_herwig[i].SetLineStyle(2)
	list_new_dy_ewk_herwig[i].SetMarkerStyle(1)
	hist_dy_ewk_herwig=list_new_dy_ewk_herwig[i].Clone("new3")
	list_new_dy_ewk_herwig[i].SetFillColor(ROOT.kViolet+6)
#	list_new_dy_ewk_herwig[i].SetFillStyle(3001)

	list_new_dy[i].Draw("E3same")
	list_new_dy_ewk_herwig[i].Draw("E3same")
	list_new_dy_ewk[i].Draw("E3same")   #CHIST

	hist_dy.Draw("CHISTsame")
	hist_dy_ewk_herwig.Draw("CHISTsame")
	hist_dy_ewk.Draw("CHISTsame")
#	list_new_dy[i].Draw("CHISTsame")
#	list_new_dy_ewk_herwig[i].Draw("CHISTsame")
#	list_new_dy_ewk[i].Draw("CHISTsame")   #CHIST
	list_new_data[i].Draw("PEsame")
	line = ROOT.TLine(xmin,1,xmax,1)
	line.SetLineStyle(3)
	line.Draw("same")
	pCMS1.Draw()
	pCMS12.Draw()
#	pCMS13.Draw()
	pCMS2.Draw()

	leg = ROOT.TLegend(0.45,0.2,0.85,0.35)
	leg.AddEntry(list_new_data[i],"Data" ,"P")
	leg.AddEntry(list_new_dy[i],"DY (%s)"%dy_model ,"L")
	leg.AddEntry(list_new_dy_ewk[i],"DY + EWK Zjj (Madgraph + Pythia8)" ,"L")
	leg.AddEntry(list_new_dy_ewk_herwig[i],"DY + EWK Zjj (Madgraph + Herwig)" ,"L")
	leg.SetFillStyle(-1)
	leg.SetBorderSize(0)
	leg.SetTextFont(42)
	leg.SetTextSize(0.025)
	leg.Draw()

	ROOT.gPad.Update()
	ROOT.gPad.RedrawAxis()
	c.SaveAs("plots/veto_eff_%s_%s.pdf"%(name,dy))
	c.SaveAs("plots/veto_eff_%s_%s.C"%(name,dy))





