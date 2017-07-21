from ROOT import TFile,TH1F
tfile = TFile.Open("../inputs/root/v25_alldata32/EWKzjj_v25_systematics.root")
ofile = TFile.Open("../inputs/root/v25_alldata32/EWKzjj_v25_systematics_range053.root","RECREATE")
new_edge_low = 0.5
new_edge_up = 3
new_hist_list = []
tfile.cd()
for h in tfile.GetListOfKeys():
    h = h.ReadObj()
    name = h.GetName()
    if name.find("atanhBDT")==-1 : continue 
#    print h.ClassName(), h.GetName()
    new_nbins = h.FindBin(new_edge_up) - h.FindBin(new_edge_low) 
 #   print new_nbins, h.FindBin(new_edge_up), h.FindBin(new_edge_low) 
   # old_up = h.GetXaxis().GetBinCenter(h.GetNbinsX())+h.GetXaxis().GetBinWidth(h.GetNbinsX())/2
    hnew = TH1F(h.GetName(),h.GetTitle(),new_nbins,new_edge_low, new_edge_up)
    for i in range(0,new_nbins) :
      hnew.SetBinContent(i+1,h.GetBinContent(i+ (h.FindBin(new_edge_low))))
      hnew.SetBinError(i+1,h.GetBinError(i+ (h.FindBin(new_edge_low))))
      #print hnew.GetBinContent(i+1), hnew.GetBinError(i+1)
    new_hist_list.append(hnew)
   # hnew.Delete()

#print len(new_hist_list)

ofile.cd()
for h in new_hist_list:
	h.Write()
