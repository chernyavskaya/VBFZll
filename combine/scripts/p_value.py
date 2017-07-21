import ROOT
from ROOT import gROOT
from ROOT import gStyle

gROOT.SetBatch(True)
gROOT.ProcessLineSync(".x /afs/cern.ch/work/n/nchernya/Hbb/setTDRStyle.C")
gROOT.ForceStyle()
gStyle.SetPadTopMargin(0.06)
gStyle.SetPadRightMargin(0.04)
gStyle.SetPadLeftMargin(0.15)

f = ROOT.TFile.Open("/mnt/t3nfs01/data01/shome/nchernya/VBFZll/combine/toys/root_mu/mlfitnumall.root")
tree = f.Get("tree_fit_sb")
hist_all = ROOT.TH1F("hist_all","",1200,-2,4)
hist = ROOT.TH1F("hist","",90,0,2)
#hist.Sumw2(True)
tree.Draw("mu>>hist","fit_status==0 && muErr<10 &&  muErr>0")
tree.Draw("mu>>hist_all","fit_status==0 && muErr<10 &&  muErr>0")
hist.Scale(1./hist.Integral())
hist_all.Scale(1./hist_all.Integral())
print hist.GetMean()
print hist.GetRMS()
hist.SetStats(False)
hist.SetLineColor(1)
hist.GetYaxis().SetTitle("Probability")
ymax=hist.GetMaximum()*1.1
hist.GetYaxis().SetRangeUser(0.,ymax)
hist.GetXaxis().SetTitle("Fitted signal strength")
c = ROOT.TCanvas("c","c",900,900)
hist.Draw("H")
func = ROOT.TF1("func","gaus",-100,100)
func.SetParameters(7.98094e-02,1.00,2.5)
func.SetLineColor(ROOT.kBlue)
func.SetLineWidth(3)
hist.Fit(func)




print func.GetChisquare(), func.GetNDF() 
#mu = 1.01696
#mu=1.02006
mu=0.995284
hist2 = hist.Clone("hist2")
hist2.SetLineColor(2)
hist2.SetLineWidth(3)
hist2.SetFillColor(ROOT.kRed-9)
hist2.SetFillStyle(1001)
hist2.GetXaxis().SetRangeUser(-100,mu)
hist2.Draw("Hsame")
#print "p_value = ", func.Integral(-100.,mu)
p_value=hist_all.Integral(hist_all.FindBin(-100.),hist_all.FindBin(mu))
print "p_value = ", p_value
line  = ROOT.TLine(1.,0.,1.,ymax)
line.SetLineColor(1)
line.SetLineStyle(2)
line.SetLineWidth(2)
line.Draw()
func2 = ROOT.TF1("func","gaus",-100,mu)
func2.FixParameter(0,func.GetParameter(0))
func2.FixParameter(1,func.GetParameter(1))
func2.FixParameter(2,func.GetParameter(2))
func2.SetLineColor(2)
func2.SetLineWidth(3)
func2.SetFillColor(ROOT.kRed-9)
func2.SetFillStyle(1001)
#func2.Draw("LFsame")
line2  = ROOT.TLine(mu,0.,mu,ymax)
line2.SetLineColor(ROOT.kRed)
line2.SetLineStyle(7)
line2.SetLineWidth(3)
line2.Draw()
pave = ROOT.TPaveText(0.7,0.8,0.9,0.9,"NDC");
pave.AddText("#mu_{obs} = %.2f" % mu)
pave.AddText("p-value = %.3f" % p_value);
pave.SetFillColor(0)
pave.SetTextFont(62)
pave.SetTextColor(2)
pave.SetTextSize(0.04)
pave.SetBorderSize(0)
pave.Draw();
pave2 = ROOT.TPaveText(0.2,0.85,0.4,0.9,"NDC");
#pave2.AddText("Combined");
#pave2.AddText("Dielectron");
pave2.AddText("Dimuon");
pave2.SetFillColor(0)
pave2.SetTextFont(42)
pave2.SetTextColor(ROOT.kBlue)
pave2.SetTextSize(0.04)
pave2.SetBorderSize(0)
pave2.Draw();
c.SaveAs("plot_pvalue_big_mu.pdf")
c.SaveAs("plot_pvalue_big_mu.png")
