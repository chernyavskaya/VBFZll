from ROOT import TFile,TH1F
tfile = TFile.Open("../inputs/root/v25_alldata5/EWKzjj_v25_systematics.root")
ofile = TFile.Open("../inputs/root/v25_alldata5/EWKzjj_v25_systematics_qgl.root","RECREATE")
new_hist_list = []
tfile.cd()
qgl_norm = {}
for h in tfile.GetListOfKeys():
    h = h.ReadObj()
    name = h.GetName()
    if name.find("atanhBDT")==-1 : continue 
    hnew = h.Clone()
    if (name.find('_QGL') == -1) and (name.find('SingleElectron') == -1) and (name.find('SingleMuon') ==-1 ) :
       current_sample = ''
       if name.find('_el')!=-1:
          if (name.find('_CMS')!=-1) : current_sample = name[name.find('el_'):name.find('_CMS')]
          else : current_sample = name[name.find('el_'):]
          qgl_norm["el_EWK_LLJJ"]=0.938595977
          qgl_norm["el_interference"]=0.938595977
          qgl_norm["el_TT"]=1.05998369
          qgl_norm["el_WW"]=0.981059169
          qgl_norm["el_WZ"]=0.956107275
          qgl_norm["el_ZZ"]=0.970928133
          qgl_norm["el_ST_tW_antitop"]=0.999659674
          qgl_norm["el_ST_tW_top"]=0.980208626
          qgl_norm["el_ST_s-channel"]=0.99449528
          qgl_norm["el_ST_t-channel_top_4f_inclusiveDecays"]=0.998267176
          qgl_norm["el_ST_t-channel_antitop_4f_inclusiveDecays"]=0.856539545
          qgl_norm["el_DYJetstoLL_amc_0J"]=0.996136594
          qgl_norm["el_DYJetstoLL_amc_1J"]=0.956949008
          qgl_norm["el_DYJetstoLL_amc_2J"]=0.952277759
          qgl_norm["el_EWKinterference"]=0.983692
       if name.find('_mu')!=-1:
          if (name.find('_CMS')!=-1) : current_sample = name[name.find('mu_'):name.find('_CMS')]
          else : current_sample = name[name.find('mu_'):]
          qgl_norm["mu_EWK_LLJJ"]=0.939774091
          qgl_norm["mu_interference"]=0.939774091
          qgl_norm["mu_TT"]=1.069615388
          qgl_norm["mu_WW"]=0.927930632
          qgl_norm["mu_WZ"]=0.967820083
          qgl_norm["mu_ZZ"]=0.964110717
          qgl_norm["mu_ST_tW_antitop"]=1.012466766
          qgl_norm["mu_ST_tW_top"]=0.990673372
          qgl_norm["mu_ST_s-channel"]=0.911006075
          qgl_norm["mu_ST_t-channel_top_4f_inclusiveDecays"]=0.986731357
          qgl_norm["mu_ST_t-channel_antitop_4f_inclusiveDecays"]=0.994925429
          qgl_norm["mu_DYJetstoLL_amc_0J"]=1.002031915
          qgl_norm["mu_DYJetstoLL_amc_1J"]=0.966710372
          qgl_norm["mu_DYJetstoLL_amc_2J"]=0.954117783
          qgl_norm["mu_EWKinterference"]=0.983032

  #     print name, current_sample, qgl_norm[current_sample]
#       print name, current_sample
       hnew.Scale(qgl_norm[current_sample]) 
    new_hist_list.append(hnew)

ofile.cd()
for h in new_hist_list:
	h.Write()
