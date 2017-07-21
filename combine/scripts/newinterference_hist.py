from ROOT import TFile,TH1F
tfile = TFile.Open("../inputs/root/v25_alldata5/EWKzjj_v25_systematics_qgl_systFull.root")
ofile = TFile.Open("../inputs/root/v25_alldata5/EWKzjj_v25_systematics_qgl_systFull_interNew.root","RECREATE")
new_hist_list = []
tfile.cd()
qgl_norm = {}
channels=['mu','el']

for channel in channels:
       h_nominal = tfile.Get("atanhBDT_%s_EWKinterference"%(channel))
       h_nominal_old = tfile.Get("atanhBDT_%s_interference"%(channel))
       h_Up_old = tfile.Get("atanhBDT_%s_interference_CMS_ewkzjj_int_shapeUp"%(channel))
       hnew_interUp = h_nominal.Clone("atanhBDT_%s_EWKinterference_CMS_ewkzjj_int_shapeUp"%(channel))
       hnew_interNom = h_nominal.Clone("atanhBDT_%s_EWKinterference"%(channel))

       for i in range(0,h_nominal.GetNbinsX()):
         if h_Up_old.GetBinContent(i+1)!=0 :
           hnew_interNom.SetBinContent(i+1, hnew_interUp.GetBinContent(i+1)*h_nominal_old.GetBinContent(i+1)/h_Up_old.GetBinContent(i+1))
           hnew_interNom.SetBinError(i+1, hnew_interUp.GetBinError(i+1)*h_nominal_old.GetBinContent(i+1)/h_Up_old.GetBinContent(i+1))
         else :     
           hnew_interNom.SetBinContent(i+1, hnew_interUp.GetBinContent(i+1))
           hnew_interNom.SetBinError(i+1, hnew_interUp.GetBinError(i+1))

       hnew_interDown = hnew_interNom.Clone("atanhBDT_%s_EWKinterference_CMS_ewkzjj_int_shapeDown"%(channel))
       hnew_interDown.Multiply(hnew_interNom)
       hnew_interDown.Divide(hnew_interUp)
       new_hist_list.append(hnew_interNom)
       new_hist_list.append(hnew_interDown)
       new_hist_list.append(hnew_interUp)


for h in tfile.GetListOfKeys():
    h = h.ReadObj()
    name = h.GetName()
    hnew = h.Clone()
    if (name.find('EWKinterference') == -1) and (name.find('interference') ==-1 ) :
       new_hist_list.append(hnew)
ofile.cd()
for h in new_hist_list:
	h.Write()
