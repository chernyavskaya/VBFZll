import ROOT
import sys

ROOT.gStyle.SetOptStat(0)


ifile = ROOT.TFile(sys.argv[1])
channel = sys.argv[2]
dy = 'amc'
datafile = ROOT.TFile("../shapes/ewkZjj_13TeV_shapes_%s_DY_%s_atanh_v25alldata5_qgl05.root" %(channel,dy))
data_hist = datafile.Get("atanhBDT_%s_data_obs" % channel)
Top = ifile.Get("shapes_fit_s/%s/Top" % channel)
EWKZ = ifile.Get("shapes_fit_s/%s/EWKZ" % channel)
VV = ifile.Get("shapes_fit_s/%s/VV" % channel)
DY = ifile.Get("shapes_fit_s/%s/DY" % channel)
interference = ifile.Get("shapes_fit_s/%s/EWKinterference" % channel)

#data = ROOT.TGraphErrors()
#for ibin in range(1,data_hist.GetNbinsX()+1):
 #  data.SetPoint(ibin, ibin-0.5, data_hist.GetBinContent(ibin))
 #  data.SetPointError(ibin,0.5, data_hist.GetBinError(ibin))
nbins=DY.GetNbinsX()
data = ROOT.TH1F("data_hist","data_hist",nbins,DY.GetBinLowEdge(1),DY.GetBinCenter(nbins)+DY.GetXaxis().GetBinWidth(nbins)/2)
for ibin in range(1,data_hist.GetNbinsX()+1):
   data.SetBinContent(ibin, data_hist.GetBinContent(ibin))
   data.SetBinError(ibin,data_hist.GetBinError(ibin))

data.SetMarkerStyle(20)
EWKZ.SetLineColor(ROOT.kRed-4)
VV.SetLineColor(ROOT.kSpring+5)
DY.SetLineColor(ROOT.kOrange-2)
Top.SetLineColor(ROOT.kBlue-4)
interference.SetLineColor(ROOT.kMagenta-9)
EWKZ.SetFillStyle(1001)
VV.SetFillStyle(1001)
DY.SetFillStyle(1001)
Top.SetFillStyle(1001)
interference.SetFillStyle(1001)
EWKZ.SetFillColor(ROOT.kRed-4)
VV.SetFillColor(ROOT.kSpring+5)
DY.SetFillColor(ROOT.kOrange-2)
Top.SetFillColor(ROOT.kBlue-4)
interference.SetFillColor(ROOT.kMagenta-9)



stack = ROOT.THStack("stack","stack")
stack.Add(VV)
stack.Add(Top)
stack.Add(DY)
stack.Add(interference)
stack.Add(EWKZ)

canv = ROOT.TCanvas("canv","canv",600,600)
canv.cd()
canv.SetBottomMargin(0.3)
xmin=0
xmax=21
frame = ROOT.TH1F("frame","",1,xmin,xmax)
canv.SetLogy()
frame.SetStats(0)
frame.GetXaxis().SetLabelSize(0)
frame.GetXaxis().SetTitleOffset(1.1);
frame.GetYaxis().SetTitleOffset(1.3);
frame.GetYaxis().SetTitle("Events")
#frame.GetXaxis().SetTitle("atanhBDT")
frame.GetYaxis().SetLabelSize(0.04)
frame.GetYaxis().SetRangeUser(10E-01,stack.GetStack().Last().GetMaximum()*100)
#frame.GetYaxis().SetRangeUser(10E-03,200)
frame.Draw()

stack.Draw("samehist")
#print "Chi2/NDOF = ",data.Chi2Test( stack.GetHistogram() , "UWCHI2/NDF")
#print "Total bkg: ",stack.GetHistogram().GetIntegral()
#stack.Draw("ep same")
stack.SetTitle(channel)
#stack.GetXaxis().SetTitle("Mjj [GeV]")
data.Draw("Psame")
#print "Total data: ",data.Integral()
#bkg.Draw("hist")
#bkg.Draw("ep same")
#stackvv.Draw("hist same")
#stackvv.Draw("ep same")
#EWKZ.Draw("HISTsame")

leg = ROOT.TLegend(0.7,0.7,0.9,0.9)
leg.SetFillStyle(-1)
leg.SetBorderSize(0)
leg.AddEntry(VV,"VV","F")
leg.AddEntry(Top,"Top","F")
leg.AddEntry(DY,"DY","F")
leg.AddEntry(interference,"interference","F")
leg.AddEntry(EWKZ,"EWKZ","F")
leg.AddEntry(data,"Data","EPL")
leg.Draw("same")

pad2 = ROOT.TPad("pad2", "pad2", 0., 0., 1., 1.)
pad2.SetTopMargin(0.73)
#	pad2.SetBottomMargin(0.)
#	pad2.SetRightMargin(0.)
pad2.SetFillColor(0)
pad2.SetFillStyle(0)
pad2.Draw()
pad2.cd()
frame2 = ROOT.TH1F("frame2","",1,xmin,xmax)
frame2.SetMinimum(-0.5)	
frame2.SetMaximum(0.5) 
frame2.GetYaxis().SetLabelSize(0.02)
frame2.GetXaxis().SetLabelSize(0.04)
frame2.GetYaxis().SetTitleSize(0.04)
frame2.GetXaxis().SetTitle("Bin number in BDT'")
frame2.SetStats(0)
frame2.GetYaxis().SetTitle("Data/MC - 1")	
frame2.Draw()

data2=data.Clone("data2")
mc_total=stack.GetStack().Last().Clone("mc_total")
data2.Add(mc_total,-1)
data2.Divide(mc_total)

mc_total_uncUp=stack.GetStack().Last().Clone("mc_total")
mc_total_uncDown=stack.GetStack().Last().Clone("mc_total")
for i in range(0,mc_total_uncUp.GetNbinsX()): 
	e=0.
	if mc_total.GetBinContent(i+1) != 0:
		e = mc_total.GetBinError(i+1)/mc_total.GetBinContent(i+1)
	mc_total_uncUp.SetBinContent(i+1,e)
	mc_total_uncDown.SetBinContent(i+1,-e)
mc_total_uncUp.SetLineColor(ROOT.kBlack)
mc_total_uncUp.SetLineWidth(1)
mc_total_uncUp.SetFillColor(ROOT.kBlack)
mc_total_uncUp.SetFillStyle(3004)
mc_total_uncDown.SetLineColor(ROOT.kBlack)
mc_total_uncDown.SetLineWidth(1)
mc_total_uncDown.SetFillColor(ROOT.kBlack)
mc_total_uncDown.SetFillStyle(3004)

mc_total_uncUp.Draw("HIST SAME")
mc_total_uncDown.Draw("HIST SAME")


line = ROOT.TLine(xmin,0,xmax,0)
line.SetLineStyle(3)
line.Draw("same")

data2.Draw("PEsame")

#print "Data yield: ",data.Integral()
#print "MC yield: ",stack.GetHistogram().Integral()

#raw_input()
#canv.SetLogy(True)
canv.SaveAs("%s_%s_postfit_alldata5.pdf" %(channel,dy))
canv.SaveAs("%s_%s_postfit_alldata5.png" %(channel,dy))
canv.SaveAs("%s_%s_postfit_alldata5.C" %(channel,dy))
ifile.Close()
datafile.Close()
