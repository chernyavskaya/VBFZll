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

f = ROOT.TFile.Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat//store/user/nchernya/VBFZll/systematics/v25/interference_%s_v25_systematics_alldata.root"%channel)

hist_names = ['atanhBDT_mu_interference','atanhBDT_mu_interference_CMS_ewkzjj_int_shapeUp','atanhBDT_mu_interference_CMS_ewkzjj_int_shapeDown']
#hist_names = ['atanhBDT_el_interference','atanhBDT_el_interference_CMS_ewkzjj_int_shapeUp','atanhBDT_el_interference_CMS_ewkzjj_int_shapeDown']
list = []

for i,name in enumerate(hist_names):
	hist = f.Get(hist_names[i])
	hist.SetLineWidth(2)
	list.append(hist)
list[0].SetLineColor(ROOT.kBlue)
list[1].SetLineColor(ROOT.kRed)
list[2].SetLineColor(ROOT.kGreen)


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


def func_norm(x,par) : 
	a0 = par[0]
	a1 = par[1]
	a2 = par[2]
	a3 = par[3]
	a4 = par[4]
	a5 = par[5]
	a6 = par[6]
	a7 = par[7]
	y = a0 + a1*x[0] + a2*pow(x[0],2)+ a3*pow(x[0],3)+ a4*pow(x[0],4)+ a5*pow(x[0],5)+ a6*pow(x[0],6)+ a7*pow(x[0],7)
	if y<0 : y=0
	return y 

def func_ratio(x,par) : 
	a0 = par[0]
	a1 = par[1]
	a2 = par[2]
	a3 = par[3]
	a4 = par[4]
	a5 = par[5]
	a6 = par[6]
	a7 = par[7]
	a02 = par[8]
	a12 = par[9]
	a22 = par[10]
	a32 = par[11]
	a42 = par[12]
	a52 = par[13]
	a62 = par[14]
	a72 = par[15]
	y1 = a0 + a1*x[0] + a2*pow(x[0],2)+ a3*pow(x[0],3)+ a4*pow(x[0],4)+ a5*pow(x[0],5)+ a6*pow(x[0],6)+ a7*pow(x[0],7)
	y2 = a02 + a12*x[0] + a22*pow(x[0],2)+ a32*pow(x[0],3)+ a42*pow(x[0],4)+ a52*pow(x[0],5)+ a62*pow(x[0],6)+ a72*pow(x[0],7)
	if y2<0 : y2=0
	if y1<0 : y1=0
	if y2!=0 : y = y1*y1/y2
	else : y=0
	return y
		


interference_func = ROOT.TF1("fnew",func_norm,5,9,8)
interference_func.FixParameter(0,-3236.73471)
interference_func.FixParameter(1,3158.70903)
interference_func.FixParameter(2,-1314.92713)
interference_func.FixParameter(3,302.849052)
interference_func.FixParameter(4,-41.69131)
interference_func.FixParameter(5,3.43119691)
interference_func.FixParameter(6,-0.15633696)
interference_func.FixParameter(7,0.00304252817)

interference_func_fixed = ROOT.TF1("fold",func_norm,5,9,8)
interference_func_fixed.FixParameter(0,-2829.38504)
interference_func_fixed.FixParameter(1,2723.93909)
interference_func_fixed.FixParameter(2,-1118.42065)
interference_func_fixed.FixParameter(3,254.03378)
interference_func_fixed.FixParameter(4,-34.4852553)
interference_func_fixed.FixParameter(5,2.79853847)
interference_func_fixed.FixParameter(6,-0.125730001)
interference_func_fixed.FixParameter(7,0.00241279903)

#interference_down = ROOT.TF1("interference_down","fnew*fnew/fold",5,9)
interference_down = ROOT.TF1("interference_down",func_ratio,5,9,16)
interference_down.FixParameter(0,-3236.73471)
interference_down.FixParameter(1,3158.70903)
interference_down.FixParameter(2,-1314.92713)
interference_down.FixParameter(3,302.849052)
interference_down.FixParameter(4,-41.69131)
interference_down.FixParameter(5,3.43119691)
interference_down.FixParameter(6,-0.15633696)
interference_down.FixParameter(7,0.00304252817)
interference_down.FixParameter(8,-2829.38504)
interference_down.FixParameter(9,2723.93909)
interference_down.FixParameter(10,-1118.42065)
interference_down.FixParameter(11,254.03378)
interference_down.FixParameter(12,-34.4852553)
interference_down.FixParameter(13,2.79853847)
interference_down.FixParameter(14,-0.125730001)
interference_down.FixParameter(15,0.00241279903)


if (1>0):
	c = ROOT.TCanvas("c","c",900,900)
	c.cd()
	c.SetBottomMargin(0.15)
#	xmin = list[0].GetXaxis().GetXmin()
#	xmax = list[0].GetXaxis().GetXmax()
	frame = ROOT.TH1F("frame","",1,5,8)
	frame.SetStats(0)
#	frame.GetXaxis().SetTitleOffset(0.91)
	frame.GetXaxis().SetTitle("log(Mqq)")
	frame.GetYaxis().SetTitle("Interference/Signal")
	frame.GetYaxis().SetLabelSize(0.04)
#	frame.GetYaxis().SetRangeUser(0.,list[1].GetMaximum()*1.2)
	frame.GetYaxis().SetRangeUser(0.,.5)
	frame.Draw()
	interference_func.Draw("SAME")
	interference_func_fixed.SetLineColor(ROOT.kBlue)
	interference_func_fixed.SetLineStyle(3)
	interference_func_fixed.Draw("SAME")
	interference_down.SetLineColor(ROOT.kGreen)
	interference_down.SetLineStyle(2)
	interference_down.Draw("SAME")
#	list[0].Draw("HISTsameE")
#	list[1].Draw("HISTsameE")
#	list[2].Draw("HISTsameE")
	pCMS1.Draw()
	pCMS12.Draw()
	pCMS2.Draw()

	leg = ROOT.TLegend(0.7,0.75,0.9,0.9)
#	leg.AddEntry(list[0],"Nominal" ,"P")
#	leg.AddEntry(list[1],"Up variation" ,"L")
#	leg.AddEntry(list[2],"Down variation" ,"L")
	leg.AddEntry(interference_func,"Nominal" ,"L")
	leg.AddEntry(interference_func_fixed,"Up variation" ,"L")
	leg.AddEntry(interference_down,"Down variation" ,"L")
	leg.SetFillStyle(-1)
	leg.SetBorderSize(0)
	leg.SetTextFont(42)
	leg.SetTextSize(0.03)
	leg.Draw()

	ROOT.gPad.Update()
	ROOT.gPad.RedrawAxis()
	c.SaveAs("plots/interference_shape_unc_%s_2.pdf"%channel)





