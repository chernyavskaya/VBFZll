import ROOT
from ROOT import gROOT
from ROOT import gStyle

def getPt(item):
	return item[0]
def getEta(item):
	return item[1]


gROOT.SetBatch(True)
gROOT.ProcessLineSync(".x /afs/cern.ch/work/n/nchernya/Hbb/setTDRStyle.C")
gROOT.ForceStyle()
gStyle.SetPadTopMargin(0.06)
gStyle.SetPadRightMargin(0.04)
gStyle.SetPadLeftMargin(0.15)

#famc = ROOT.TFile.Open("/afs/cern.ch/work/n/nchernya/VBFZll/plotter/output_root/DYJetstoLL_amc_mu_QCDScalenom_JESnom_v24_ewk_mucorr_bdt.root")
#fmdg = ROOT.TFile.Open("/afs/cern.ch/work/n/nchernya/VBFZll/plotter/output_root/DYJetstoLL_HT_mu_QCDScalenom_JESnom_v24_ewk_mucorr_bdt.root")
famc = ROOT.TFile.Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat//store/user/nchernya/VBFZll/plotterOutput/v25/DYJetstoLL_amc_el_QCDScalenom_JESnom_v25_bdt_alldata.root")
fmdg = ROOT.TFile.Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat//store/user/nchernya/VBFZll/plotterOutput/v25/DYJetstoLL_HT_el_QCDScalenom_JESnom_v25_bdt_noMqq.root")
famc2 = ROOT.TFile.Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat//store/user/nchernya/VBFZll/plotterOutput/v25/DYJetstoLL_amc_mu_QCDScalenom_JESnom_v25_bdt_alldata.root")
fmdg2 = ROOT.TFile.Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat//store/user/nchernya/VBFZll/plotterOutput/v25/DYJetstoLL_HT_mu_QCDScalenom_JESnom_v25_bdt_noMqq.root")

hist_names = ['hHT','hEtaQQ',"hMqq","hMqq_log","hlheHT_log","hPhiQQ","hJet1q_pt", "hJet2q_pt", "hJets12_pt","hJets12_pt_log"]
list_amc = []
list_amc_ratio = []
list_mdg = []
#for i,name in enumerate(hist_names):
#	hist_amc = famc.Get(hist_names[i])
#	hist_mdg = fmdg.Get(hist_names[i])
#	hist_amc.SetLineColor(ROOT.kRed)
#	hist_mdg.SetLineColor(ROOT.kBlue)
#	hist_amc.SetLineWidth(2)
#	hist_mdg.SetLineWidth(2)
#	hist_mdg.SetLineStyle(2)
#	list_amc.append(hist_amc)
#	list_amc_ratio.append(hist_amc)
#	list_mdg.append(hist_mdg)

for i,name in enumerate(hist_names):
	hist_amc = famc.Get(hist_names[i])
	hist_mdg = fmdg.Get(hist_names[i])
	hist_amc2 = famc2.Get(hist_names[i])
	hist_mdg2 = fmdg2.Get(hist_names[i])
	if i==3 :
	#	new_nbins = hist_amc.FindBin(9) - hist_amc.FindBin(5.3) 
		new_nbins=97
		hnew_amc = ROOT.TH1F(hist_amc.GetName(),hist_amc.GetTitle(),new_nbins,5.3,15)
		hnew_mdg = ROOT.TH1F(hist_mdg.GetName(),hist_mdg.GetTitle(),new_nbins,5.3, 15)
		hnew_amc2 = ROOT.TH1F(hist_amc2.GetName(),hist_amc2.GetTitle(),new_nbins,5.3,15)
		hnew_mdg2 = ROOT.TH1F(hist_mdg2.GetName(),hist_mdg2.GetTitle(),new_nbins,5.3,15)
		for i in range(0,new_nbins) :
			hnew_amc.SetBinContent(i+1,hist_amc.GetBinContent(i+ (hist_amc.FindBin(5.3))))
			hnew_amc.SetBinError(i+1,hist_amc.GetBinError(i+ (hist_amc.FindBin(5.3)) ))	
			hnew_amc2.SetBinContent(i+1,hist_amc2.GetBinContent(i+ (hist_amc2.FindBin(5.3))))
			hnew_amc2.SetBinError(i+1,hist_amc2.GetBinError(i+ (hist_amc2.FindBin(5.3)) ))	
			hnew_mdg.SetBinContent(i+1,hist_mdg.GetBinContent(i+ (hist_mdg.FindBin(5.3))))
			hnew_mdg.SetBinError(i+1,hist_mdg.GetBinError(i+ (hist_mdg.FindBin(5.3)) ))	
			hnew_mdg2.SetBinContent(i+1,hist_mdg2.GetBinContent(i+ (hist_mdg2.FindBin(5.3))))
			hnew_mdg2.SetBinError(i+1,hist_mdg2.GetBinError(i+ (hist_mdg2.FindBin(5.3)) ))
		hnew_amc.GetXaxis().SetTitle(hist_amc.GetXaxis().GetTitle())	
		hist_amc=hnew_amc	
		hist_mdg=hnew_mdg
		hist_amc2=hnew_amc2	
		hist_mdg2=hnew_mdg2
#	hist_amc.SetLineColor(ROOT.kRed)
#	hist_mdg.SetLineColor(ROOT.kBlue)
#	hist_amc.SetLineWidth(2)
#	hist_mdg.SetLineWidth(2)
#	hist_mdg.SetLineStyle(2)
	hist_amc2.SetLineColor(ROOT.kRed)
	hist_mdg2.SetLineColor(ROOT.kBlue)
	hist_amc2.SetLineWidth(2)
	hist_mdg2.SetLineWidth(2)
	hist_mdg2.SetLineStyle(2)
#	hist_amc.Add(hist_amc2)
#	hist_mdg.Add(hist_mdg2)
	list_amc.append(hist_amc2)
	list_amc_ratio.append(hist_amc2)
	list_mdg.append(hist_mdg2)


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
	list_amc[num].Scale(1./list_amc[num].Integral())
	list_mdg[num].Scale(1./list_mdg[num].Integral())
	hist_ratio = list_amc_ratio[num].Clone("ratio")
	hist_ratio.Divide(list_mdg[num])

	hist_ratio_up = hist_ratio.Clone("up")	
	hist_ratio_down = hist_ratio.Clone("down")
	for i in range(0,hist_ratio.GetNbinsX()) :
		ey=hist_ratio.GetBinError(i+1)
		bin=hist_ratio.GetBinContent(i+1)
		hist_ratio_up.SetBinContent(i+1,bin+ey)	
		hist_ratio_down.SetBinContent(i+1,bin-ey)

	func_res_up = ROOT.TF1("funcup","pol6",5.3,9.)
	func_res_up.SetLineWidth(1)
	func_res_up.SetLineColor(ROOT.kRed-7) 
	func_res_up.SetFillColor(ROOT.kRed-7)  	
	func_res_up.SetFillStyle(3004)
	func_res_down = ROOT.TF1("funcdown","pol6",5.3,9.)
	func_res_down.SetLineWidth(1)
	func_res_down.SetLineColor(ROOT.kRed-7) 
	func_res_down.SetFillColor(10)  
	func_res_down.SetFillStyle(1001)
	if (num==3) :
		hist_ratio_up.Fit(func_res_up,"R","N")
		hist_ratio_down.Fit(func_res_down,"R","N")
		for i in range(0,7):
			print 'func_Mqq_up->FixParameter(%i,%f);'% (i,func_res_up.GetParameter(i))
			print 'func_Mqq_down->FixParameter(%i,%f);'% (i,func_res_down.GetParameter(i))
	c = ROOT.TCanvas("c","c",900,900)
	c.cd()
	c.SetBottomMargin(0.3)
	xmin = list_amc[num].GetXaxis().GetXmin()
	xmax = list_amc[num].GetXaxis().GetXmax()
	if num==3 : 
		xmin = 5
		xmax= 9
	frame = ROOT.TH1F("frame","",1,xmin,xmax)
#	if num==3 : 
#	frame.SetMinimum(0.00005)
#	frame.GetYaxis().SetRangeUser(1e-07,list_amc[num].GetMaximum()*1.2)
#	c.SetLogy()
	frame.SetStats(0)
	frame.GetXaxis().SetLabelSize(0)
	frame.GetXaxis().SetTitleOffset(0.91);
	frame.GetYaxis().SetTitle("Events")
	frame.GetYaxis().SetLabelSize(0.04)
#	list_amc[num].GetXaxis().SetTitle(list_amc[num].GetTitle())
	frame.GetYaxis().SetRangeUser(0.,list_amc[num].GetMaximum()*1.2)
	frame.Draw()
	list_amc[num].Draw("HISTsameE")
	list_mdg[num].Draw("HISTsameE")
	pCMS1.Draw()
	pCMS12.Draw()
	pCMS2.Draw()

	leg = ROOT.TLegend(0.7,0.75,0.9,0.9)
	#leg.AddEntry(list_amc[num],"Dielectron" ,"P")
	leg.AddEntry(list_amc[num],"Dimuon" ,"P")
#	leg.AddEntry(list_amc[num],"Dileptons" ,"P")
	leg.AddEntry(list_amc[num],"AMC@NLO" ,"L")
	leg.AddEntry(list_mdg[num],"Madgraph" ,"L")
	leg.SetFillStyle(-1)
	leg.SetBorderSize(0)
	leg.SetTextFont(42)
	leg.SetTextSize(0.03)
	leg.Draw()

	pad2 = ROOT.TPad("pad2", "pad2", 0., 0., 1., 1.)
	pad2.SetTopMargin(0.73)
#	pad2.SetBottomMargin(0.)
#	pad2.SetRightMargin(0.)
	pad2.SetFillColor(0)
	pad2.SetFillStyle(0)
	pad2.Draw()
	pad2.cd()
	frame2 = ROOT.TH1F("frame2","",1,xmin,xmax)
	frame2.SetMinimum(0.5)	
	frame2.SetMaximum(1.5) 
	frame2.GetYaxis().SetLabelSize(0.02)
	frame2.GetXaxis().SetLabelSize(0.04)
	frame2.GetYaxis().SetTitleSize(0.04)
	frame2.GetXaxis().SetTitle(list_amc[num].GetXaxis().GetTitle())
	frame2.SetStats(0)
	frame2.GetYaxis().SetTitle("AMC/MDG")	
	frame2.Draw()

	func_res_up.Draw("sAME")
	func_res_down.Draw("sAME")
	hist_ratio.Draw("HISTsameE")
	#func_res = ROOT.TF1("func","([0]+[1]*x+[2]*x*x+[3]*x*x*x)*TMath::Exp(-1*[4]*x)",4.4,7.4)
#	func_res = ROOT.TF1("func","[0]+[1]*x+[2]*x*x+[3]*x*x*x",0.,8)
#	func_res = ROOT.TF1("func","pol5",4.3,7.9)   #for muons for ptSum corr 
#	func_res.SetParameters(-460,408,-141,24,-2,0.2)  #for muons for ptSum corr
#	func_res = ROOT.TF1("func","pol5",4.3,8)   #for electrons for ptSum corr
#	func_res.SetParameters(1,-1,0.5) #for electrons for ptSum corr
#	func_res = ROOT.TF1("func","pol7",0,5.7) # for electrons EtaQQ corr 
#	func_res = ROOT.TF1("func","pol7",0,8.3) # for muons EtaQQ corr 
#	func_res = ROOT.TF1("func","pol5",0,8.5) # for electrons EtaQQ corr 
#	func_res = ROOT.TF1("func","pol7",0,8.3) # for muons EtaQQ corr 
#	func_res = ROOT.TF1("func","pol6",5.32,8.3) # for muons Mqq_log
#	func_res = ROOT.TF1("func","pol3",7.55,8.0) # for electrons Mqq_log
#	func_res_copy = ROOT.TF1("func2","pol3",7.55,8.0) # for electrons Mqq_log
#	func_res2 = ROOT.TF1("func_copy","pol4",5.2,7.55) # for electrons Mqq_log
#	func_res = ROOT.TF1("func","pol6",5.2,8.3) # for bpth leptons Mqq_log
	func_res = ROOT.TF1("func","pol6",5.3,9) # for bpth leptons Mqq_log
#	func_res.SetParameters(861,-778,292,-58,6.5,-0.4,0.01)
	func_res.SetLineColor(ROOT.kGreen-3)


#	func_res_copy.SetLineColor(ROOT.kGreen-3)
#	func_res2.SetLineColor(ROOT.kGreen-3)
#	if (num==len(hist_names)-1) :
	if (num==3) :
		hist_ratio.Fit(func_res,"R","SAMe")
#		for i in range(0,4):
#			func_res_copy.FixParameter(i,func_res.GetParameter(i))
#		hist_ratio.Fit(func_res2,"R","SAMe")
#		func_res_copy.Draw("L SAME")
		for i in range(0,7):
			print 'func_Mqq->FixParameter(%i,%f);'% (i,func_res.GetParameter(i))
#		print hist_names[num],func_res.GetChisquare(), func_res.GetNDF(), func_res.Eval(7.7113)
		print hist_names[num],func_res.GetChisquare(), func_res.GetNDF(), func_res.GetChisquare()/func_res.GetNDF(), func_res.GetX(1.,5.1,5.3),func_res.GetX(1.,7.3,9)
	line = ROOT.TLine(xmin,1,xmax,1)
	line.SetLineStyle(3)
	line.Draw("same")
	ROOT.gPad.Update()
	ROOT.gPad.RedrawAxis()
#	c.SaveAs("plots_amc_mdg_corr/plot_amc_mdg_corr_mu_Ptsumcorr_ForEta_%s_after.pdf" %list_amc[num].GetName() )
	c.SaveAs("plots_amc_mdg_corr/v25/plot_newest4_amc_mdg_corr_mu_%s.pdf" %list_amc[num].GetName() )





