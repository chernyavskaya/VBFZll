import ROOT
import array
import os
import json

# A class to apply SF's tabulated in json files
class LeptonSF:
    def __init__(self, lep_json, lep_name, lep_binning, extrapolateFromClosestBin=True) :
        if not os.path.isfile(lep_json):
            self.valid = False
            if lep_json!="":
                print "[LeptonSF]: Warning: ", lep_json, " is not a valid file. Return."
            else:
                print "[LeptonSF]: No file has been specified. Return."
        else:
            self.init(lep_json, lep_name, lep_binning, extrapolateFromClosestBin)

    def init(self, lep_json, lep_name, lep_binning, extrapolateFromClosestBin) :
        f = open(lep_json, 'r')             
        print '[LeptonSF]: Initialize with the following parameters:'
        print '\tfile:',lep_json
        print '\titem:', lep_name
        print '\tbinning:', lep_binning
        results = json.load(f)
        if lep_name not in results.keys():
            self.valid = False
            print "[LeptonSF]: Warning: ", lep_name , " is not a valid item. Return."
            return False
        self.res = results[lep_name]
        self.lep_name = lep_name
        self.lep_binning = lep_binning
        self.valid = True
        self.extrapolateFromClosestBin = extrapolateFromClosestBin
        f.close()

    def get_1D(self, pt):
        if not self.valid:
            return [1.0, 0.0]

        stripForEta = 5
        if self.lep_binning not in self.res.keys():
            return [1.0, 0.0]

        # if no bin is found, search for closest one, and double the uncertainty
        closestPtBin = ""
        closestPt = 9999.

        ptFound = False

        for ptKey, result in sorted(self.res[self.lep_binning].iteritems()) :
            #print 'ptKey is', ptKey
            ptL = float(((ptKey[7:]).rstrip(']').split(',')[0]))
            ptH = float(((ptKey[7:]).rstrip(']').split(',')[1]))

            #print 'ptL is', ptL
            #print 'ptH is', ptH

            if abs(ptL-pt)<closestPt or abs(ptH-pt)<closestPt and not ptFound:
                closestPt = min(abs(ptL-pt), abs(ptH-pt))
                closestPtBin = ptKey

                if (pt>ptL and pt<ptH):
                    closestPtBin = ptKey
                    ptFound = True

                if ptFound:
                    return [result["value"], result["error"]]

        if self.extrapolateFromClosestBin and not (closestPtBin==""):
            return [self.res[self.lep_binning][closestPtBin]["value"],2*self.res[self.lep_binning][closestPtBin]["error"]]
        else:
            return [1.0, 0.0]
                    

    def get_2D(self, pt, eta):
        if not self.valid:
            return [1.0, 0.0]        

        stripForEta = 5
        if self.lep_binning not in self.res.keys():
            return [1.0, 0.0]

        if "abseta" in self.lep_binning:
            eta = abs(eta)
            stripForEta = 8

        # if no bin is found, search for closest one, and double the uncertainty
        closestEtaBin = ""
        closestPtBin = ""
        closestEta = 9999.
        closestPt = 9999.

        etaFound = False
        for etaKey, values in sorted(self.res[self.lep_binning].iteritems()) :
            etaL = float(((etaKey[stripForEta:]).rstrip(']').split(',')[0]))
            etaH = float(((etaKey[stripForEta:]).rstrip(']').split(',')[1]))

            ptFound = False

            if abs(etaL-eta)<closestEta or abs(etaH-eta)<closestEta and not etaFound:
                closestEta = min(abs(etaL-eta), abs(etaH-eta))
                closestEtaBin = etaKey

            if (eta>etaL and eta<etaH):
                closestEtaBin = etaKey
                etaFound = True                

            for ptKey, result in sorted(values.iteritems()) :
                ptL = float(((ptKey[4:]).rstrip(']').split(',')[0]))
                ptH = float(((ptKey[4:]).rstrip(']').split(',')[1]))                

                if abs(ptL-pt)<closestPt or abs(ptH-pt)<closestPt and not ptFound:
                    closestPt = min(abs(ptL-pt), abs(ptH-pt))
                    closestPtBin = ptKey

                if (pt>ptL and pt<ptH):
                    closestPtBin = ptKey
                    ptFound = True

                if etaFound and ptFound:
                    return [result["value"], result["error"]]

        if self.extrapolateFromClosestBin and not (closestPtBin=="" or closestEtaBin==""):
            return [self.res[self.lep_binning][closestEtaBin][closestPtBin]["value"], 
                    2*self.res[self.lep_binning][closestEtaBin][closestPtBin]["error"]] 
        else:
            return [1.0, 0.0]


    def triggerMap(self):
        if not self.valid:
		      return -1       
        stripForEta = 5
        if self.lep_binning not in self.res.keys():
				return -1

        if "abseta" in self.lep_binning:
            stripForEta = 8

        etaBins = []
        ptBins = []

        # if no bin is found, search for closest one, and double the uncertainty
        closestEtaBin = ""
        closestPtBin = ""
        closestEta = 9999.
        closestPt = 9999.

        etaFound = False
        for etaKey, values in sorted(self.res[self.lep_binning].iteritems()) :
            etaL = float(((etaKey[stripForEta:]).rstrip(']').split(',')[0]))
            etaH = float(((etaKey[stripForEta:]).rstrip(']').split(',')[1]))
            etaBins.append(etaL) 
            etaBins.append(etaH)	
            for ptKey, result in sorted(values.iteritems()) :
                ptL = float(((ptKey[4:]).rstrip(']').split(',')[0]))
                ptH = float(((ptKey[4:]).rstrip(']').split(',')[1]))             
                ptBins.append(ptL)
                ptBins.append(ptH)

        etaBins = list(set(etaBins)) # get rid of duplicates
        ptBins = list(set(ptBins))
        etaBins = sorted(etaBins) 
        ptBins = sorted(ptBins)

        triggerEffMap = ROOT.TH2F("TriggerEffMap_"+self.lep_name,"TriggerEffMap_"+self.lep_name,len(ptBins)-1,ptBins[0],ptBins[len(ptBins)-1],len(etaBins)-1,etaBins[0],etaBins[len(etaBins)-1] )
       # triggerEffMap = ROOT.TH2F("tt","tt", len(ptBins)-1, array(ptBins), len(etaBins)-1, array(etaBins))
        #print array(ptBins), array(etaBins)
        
        for etaKey, values in sorted(self.res[self.lep_binning].iteritems()):
                etaL = float(((etaKey[stripForEta:]).rstrip(']').split(',')[0]))
                etaH = float(((etaKey[stripForEta:]).rstrip(']').split(',')[1]))
                binLow1 = triggerEffMap.GetYaxis().FindBin(etaL)
                binHigh1 = triggerEffMap.GetYaxis().FindBin(etaH) - 1
                for ptKey, result in sorted(values.iteritems()):
                    ptL = float(((ptKey[4:]).rstrip(']').split(',')[0]))
                    ptH = float(((ptKey[4:]).rstrip(']').split(',')[1]))
                    binLow2 = triggerEffMap.GetXaxis().FindBin(ptL)
                    binHigh2 = triggerEffMap.GetXaxis().FindBin(ptH) - 1
                    for ibin1 in range(binLow1,binHigh1+1):
                        for ibin2 in range(binLow2,binHigh2+1):
                            triggerEffMap.SetBinContent(ibin2,ibin1, result["value"])
                            triggerEffMap.SetBinError(ibin2,ibin1, result["error"])
                            print result["value"]
		
        fout = ROOT.TFile("TriggerEffMap_"+self.lep_name+".root","RECREATE")
        fout.cd()
        triggerEffMap.Write()
        fout.Close()
        return 0

	

##################################################################################################
# EXAMPLE 
#

if __name__ == "__main__":

  #  jsonpath = os.environ['CMSSW_BASE']+"/src/VHbbAnalysis/Heppy/data/leptonSF/"
    jsonpath = "/afs/cern.ch/work/n/nchernya/VBFZll/plotter/trigger_efficiens/"
    jsons = {    
     #   jsonpath+'SingleMuonTrigger_LooseMuons_beforeL2fix_Z_RunBCD_prompt80X_7p65.json' :['MuonTrigger_data_all_IsoMu22_OR_IsoTkMu22_pteta_Run2016B_beforeL2Fix', 'abseta_pt_MC'],
    #    jsonpath+'eff_Ele27_WPLoose_Eta2p1_RunBCD.json' :['HLT_Ele27_WPLoose_eta2p1_WP80_BCD', 'eta_pt_ratio'],
        #jsonpath+'SingleMuonTrigger_Z_RunBCD_prompt80X_7p65.json' :['IsoMu22_OR_IsoTkMu22_PtEtaBins_Run273158_to_274093', 'abseta_pt_DATA' ],
        #jsonpath+'MuonIso_Z_RunBCD_prompt80X_7p65.json' : ['MC_NUM_LooseRelIso_DEN_TightID_PAR_pt_spliteta_bin1', 'abseta_pt_ratio'],
        #jsonpath+'SingleMuonTrigger_Z_RunBCD_prompt80X_7p65.json' :['IsoMu22_OR_IsoTkMu22_PtEtaBins_Run274094_to_276097', 'abseta_pt_DATA' ],
        #jsonpath+'MuonTrkHIP_80X_Jul28.json' :[ 'ratio_eta', 'ratio_eta' ],
        #jsonpath+'MuonTrkHIP_80X_Jul28.json' :['ratio_vtx', 'ratio_vtx' ],
        #jsonpath+'SingleMuonTrigger_Z_RunCD_Reco76X_Feb15.json' : ['runC_IsoMu20_OR_IsoTkMu20_PtEtaBins', 'abseta_pt_ratio' ]
        #jsonpath+'SingleMuonTrigger_Z_RunCD_Reco76X_Feb15.json' : ['runD_IsoMu20_OR_IsoTkMu20_HLTv4p2_PtEtaBins', 'abseta_pt_ratio' ],
     #   jsonpath+'SingleMuonTrigger_LooseMuons_afterL2fix_Z_RunBCD_prompt80X_7p65.json' : ['MuonTrigger_data_all_IsoMu22_OR_IsoTkMu22_pteta_Run2016B_afterL2Fix', 'abseta_pt_MC' ],
     #   jsonpath+'SingleMuonTrigger_LooseMuons_beforeL2fix_Z_RunBCD_prompt80X_7p65.json' : ['MuonTrigger_data_all_IsoMu22_OR_IsoTkMu22_pteta_Run2016B_beforeL2Fix', 'abseta_pt_MC' ],
     #   jsonpath+'WP90_BCDEF_withRelIso.json' : ['electronTriggerEfficiencyHLT_Ele27_WPLoose_eta2p1_WP90_BCDEF','eta_pt_ratio'],
     #   jsonpath+'SingleMuonTrigger_BCDEF.json' : ['IsoMu24_OR_IsoTkMu24_PtEtaBins', 'pt_abseta_ratio' ],
      #  jsonpath+'SingleMuonTrigger_GH.json' : ['IsoMu24_OR_IsoTkMu24_PtEtaBins', 'pt_abseta_ratio' ],
        #jsonpath+'SingleMuonId_BCDEF.json' : ['MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta', 'pt_abseta_ratio' ],
        #jsonpath+'SingleMuonId_GH.json' : ['MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta', 'pt_abseta_ratio' ],
      #  jsonpath+'SingleMuonIso_BCDEF.json' : ['LooseISO_LooseID_pt_eta', 'pt_abseta_ratio' ],
        jsonpath+'SingleMuonIso_GH.json' : ['LooseISO_LooseID_pt_eta', 'pt_abseta_ratio' ],
    #    jsonpath+'Electron_Id_WP80.json' : ['ScaleFactor_MVAIDWP80_80x', 'eta_pt_ratio' ],
     #   jsonpath+'Electron_Id_WP90.json' : ['ScaleFactor_MVAID_80x', 'eta_pt_ratio' ],
      #  jsonpath+'Electron_tracker.json' : ['ScaleFactor_tracker_80x', 'eta_pt_ratio' ],
        }
    for j, name in jsons.iteritems():
        lepCorr = LeptonSF(j , name[0], name[1])
        lepCorr.triggerMap()
      #  weight = lepCorr.get_2D( 65 , -1.5)
      #  val = weight[0]
      #  err = weight[1]
      #  print 'SF: ',  val, ' +/- ', err
    
    
    #jsons = {
    #    'SingleMuonTrigger_Z_RunCD_Reco74X_Dec1.json' : ['runD_IsoMu20_OR_IsoTkMu20_HLTv4p3_PtEtaBins', 'abseta_pt_ratio'],
    #    'MuonIso_Z_RunCD_Reco74X_Dec1.json' : ['NUM_LooseRelIso_DEN_LooseID_PAR_pt_spliteta_bin1', 'abseta_pt_ratio'], 
    #    'MuonID_Z_RunCD_Reco74X_Dec1.json' : ['NUM_LooseID_DEN_genTracks_PAR_pt_spliteta_bin1', 'abseta_pt_ratio'] ,
    #    'CutBasedID_LooseWP.json' : ['CutBasedID_LooseWP', 'eta_pt_ratio'],
    #    'CutBasedID_TightWP.json' : ['CutBasedID_TightWP', 'eta_pt_ratio'],
    #    'SingleMuonTrigger_Z_RunCD_Reco74X_Dec1_MC.json' : ['runD_IsoMu20_OR_IsoTkMu20_HLTv4p3_PtEtaBins', 'abseta_pt_MC'],
    #    'SingleMuonTrigger_Z_RunCD_Reco74X_Dec1_MC.json' : ['runD_IsoMu20_OR_IsoTkMu20_HLTv4p3_PtEtaBins', 'abseta_pt_MC'],
    #    }
    #
    #for j, name in jsons.iteritems():
    #    lepCorr = LeptonSF(j , name[0], name[1])
    #    weight = lepCorr.get_2D( 35. , 1.0)
    #    val = weight[0]
    #    err = weight[1]
    #    print j, name[0], ': ',  val, ' +/- ', err
