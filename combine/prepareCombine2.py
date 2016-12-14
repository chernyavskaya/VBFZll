#!/usr/bin/env python
from ROOT import TFile, TH1, TH1D, TAxis
 
tag13var = "20140414_Fast_13var"

plotCuts = ["_200","_ptHard_200","_ptHard_ystar_200"]

mcGroups = {'DY_mdg': ["DYJetstoLL_HT100","DYJetstoLL_HT100_200","DYJetstoLL_HT200_400","DYJetstoLL_HT400_600","DYJetstoLL_HT600_Inf"],
				'DY_amc': ["DYJetstoLL_amc"],
            'Top': ["TT"],
            'VV': ["WW","WZ","ZZ"],
            'EWKZ': ["EWK_LLJJ"]}

def getTreeLocation():
  return "/afs/cern.ch/work/n/nchernya/VBFZll/combine/"

def addErrorToHist(h, sigma, bin = 0):
  if sigma == 0: return h;
  h.SetBinContent(bin, h.GetBinContent(bin) + h.GetBinError(bin)*sigma);
  return h;

def getPlot(sourceFile, sample, plot, ignoreNonExist = False, addError = 0., bin = 0):
  hist = sourceFile.Get(sample+"/"+plot)
  if sample=="" : hist=sourceFile.Get(plot)
  if hist and hist.GetEntries() != 0: return addErrorToHist(hist.Clone(), addError, bin)
  if not ignoreNonExist: 
    print plot + " not found!"
    exit(1)

def safeAdd(first, second,name):
  if first is None: 
    print "first doesn't exist"
    return second
  if second is None:
    print "second doesn't exist"
    return first
  merged = first.Clone(name)
  merged.Add(second)
  return merged

def merge(sourceFile, plot, histList, addError = 0., bin = 0):
  histMerged = getPlot(sourceFile, histList[0], plot, False, addError, bin)
  for hist in histList[1:]: histMerged = safeAdd(histMerged, getPlot(sourceFile, hist, plot, False, addError, bin),plot)
  return histMerged

expected = {}
expectedStr = {}
JESnorm = {}
combineFile = TFile("ewkZjj_13TeV_2.root","RECREATE")
for type in ["mu","el"]: 
 for process in ["SingleMuon","SingleElectron"] : #,"WW","WZ","ZZ","TT","EWK_LLJJ","DYJetstoLL_amc","DYJetstoLL_HT100","DYJetstoLL_HT100_200","DYJetstoLL_HT200_400","DYJetstoLL_HT400_600","DYJetstoLL_HT600_Inf"]:
  if process=="SingleMuon" and type!="mu" :
    continue
  if process=="SingleElectron" and type!="el" :
    continue
  print "prepare " + type + process
  sourceFile13var = TFile(getTreeLocation() + "ewkZjj_13TeV.root")
  for basePlot in ["BDT","atanhBDT"]:
    sourceFile = sourceFile13var
    basePlot_ = basePlot
    data = getPlot(sourceFile, "", basePlot_ +"_"+ type +"_"+ process )
    combineFile.cd()
    data.Write()
 
  for basePlot in ["BDT","atanhBDT"]:
    sourceFile = sourceFile13var
    basePlot_ = basePlot
    sourceFile.cd()
    for systematic in [""]:
      for name, mcs in mcGroups.iteritems():
        plot = basePlot_ +"_"+ type +"_"+ name
        thisGroup = merge(sourceFile, plot, mcs)
        combineFile.cd()
        thisGroup.Write(plot)

exit()
