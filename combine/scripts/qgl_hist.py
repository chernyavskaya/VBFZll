from ROOT import TFile,TH1F
tfile = TFile.Open("../inputs/root/v25_alldata5/EWKzjj_v25_systematics_qgl.root")
ofile = TFile.Open("../inputs/root/v25_alldata5/EWKzjj_v25_systematics_qgl_systFull.root","RECREATE")
names = ['EWK_LLJJ','DYJetstoLL_amc_0J','DYJetstoLL_amc_1J','DYJetstoLL_amc_2J','ZZ','WZ','WW','TT','ST_t-channel_top_4f_inclusiveDecays','ST_tW_antitop','ST_s-channel','ST_tW_top','ST_t-channel_antitop_4f_inclusiveDecays']
new_hist_list = []
tfile.cd()
qgl_norm = {}
channels=['mu','el']

for channel in channels:
    for i in range(0,len(names)):
       h_nominal = tfile.Get("atanhBDT_%s_%s"%(channel,names[i]))
       h_qglUp = (tfile.Get("atanhBDT_%s_%s_CMS_ewkzjj_QGLUp"%(channel,names[i]))).Clone("Up")
       h_qglDown = tfile.Get("atanhBDT_%s_%s_CMS_ewkzjj_QGLDown"%(channel,names[i])).Clone("Down")
       hnew_qglUp = h_nominal.Clone("atanhBDT_%s_%s_CMS_ewkzjj_QGLUp"%(channel,names[i]))
       hnew_qglUp.Multiply(h_nominal)
       hnew_qglUp.Divide(h_qglUp)
       hnew_qglDown = h_qglUp.Clone("atanhBDT_%s_%s_CMS_ewkzjj_QGLDown"%(channel,names[i]))
       new_hist_list.append(hnew_qglUp)
       new_hist_list.append(hnew_qglDown)

    #   hnew_qglDown = h_qglDown.Clone("atanhBDT_%s_%s_CMS_ewkzjj_QGLDown"%(channel,names[i]))
     #  hnew_qglDown.Add(h_nominal,hnew_qglDown,0.75,0.25)
     #  hnew_qglUp = h_nominal.Clone("atanhBDT_%s_%s_CMS_ewkzjj_QGLUp"%(channel,names[i]))
     #  hnew_qglUp.Multiply(h_nominal)
     #  hnew_qglUp.Divide(hnew_qglDown)
     #  new_hist_list.append(hnew_qglUp)
     #  new_hist_list.append(hnew_qglDown)


for h in tfile.GetListOfKeys():
    h = h.ReadObj()
    name = h.GetName()
    hnew = h.Clone()
    if (name.find('_QGL') == -1) :
       new_hist_list.append(hnew)
ofile.cd()
for h in new_hist_list:
	h.Write()
