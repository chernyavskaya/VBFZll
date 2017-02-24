import ROOT
import sys

ROOT.gStyle.SetOptStat(0)


channel = sys.argv[1]
dy = sys.argv[2]
ifile = ROOT.TFile("../shapes/ewkZjj_13TeV_shapes_%s_DY_%s_atanh_v25.root" %(channel,dy))
data = ifile.Get("atanhBDT_%s_data_obs" % channel)
Top = ifile.Get("atanhBDT_%s_Top" % channel)
EWKZ = ifile.Get("atanhBDT_%s_EWKZ" % channel)
VV = ifile.Get("atanhBDT_%s_VV" % channel)
DY = ifile.Get("atanhBDT_%s_DY" % channel)
interference = ifile.Get("atanhBDT_%s_interference" % channel)


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

canv = ROOT.TCanvas("canv","canv")
frame = ROOT.TH1F("frame","",1,0,3)
canv.SetLogy()
frame.SetStats(0)
frame.GetXaxis().SetLabelSize(0)
frame.GetXaxis().SetTitleOffset(0.91);
frame.GetYaxis().SetTitle("Events")
frame.GetXaxis().SetTitle("atanhBDT")
frame.GetYaxis().SetLabelSize(0.04)
frame.GetYaxis().SetRangeUser(10E-03,data.GetMaximum()*100)
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
leg.AddEntry(VV,"VV")
leg.AddEntry(Top,"Top")
leg.AddEntry(DY,"DY")
leg.AddEntry(interference,"interference")
leg.AddEntry(EWKZ,"EWKZ")
leg.AddEntry(data,"Data")
leg.Draw("same")

#print "Data yield: ",data.Integral()
#print "MC yield: ",stack.GetHistogram().Integral()

#raw_input()
#canv.SetLogy(True)
canv.SaveAs("%s_%s_prefit.pdf" %(channel,dy))
canv.SaveAs("%s_%s_prefit.png" %(channel,dy))
canv.SaveAs("%s_%s_prefit.C" %(channel,dy))
ifile.Close()


