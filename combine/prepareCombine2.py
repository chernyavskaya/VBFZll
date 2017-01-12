#!/usr/bin/env python
from ROOT import TFile, TH1,TH1F, TH1D, TAxis
import numpy
from math import sqrt
import os,sys 
from optparse import OptionParser,OptionGroup

tempargv = sys.argv[:]
sys.argv = []
sys.argv = tempargv


directoryName = '/afs/cern.ch/work/n/nchernya/VBFZll/combine/shapes/'
NBINS=25

def parser(mp=None): 
  if not mp: mp = OptionParser()
  mp.add_option('--dy',help='DY sample used, mdg or amc',default='amc',type='str',dest='dy_sample')
  mp.add_option('--bdttype',help='BDT or transformed BDT',default='atanh',type='str',dest='bdt_type')
  return mp
 

def getTreeLocation():
  return "/afs/cern.ch/work/n/nchernya/VBFZll/combine/"

def addErrorToHist(h, sigma, bin = 0):
  hist = h.Clone()
  if sigma == 0: return h;
  hist.SetBinContent(bin, h.GetBinContent(bin) + h.GetBinError(bin)*sigma);
  return hist;

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


mp = parser()
opts,args = mp.parse_args()
dy_option = opts.dy_sample
bdt_option = opts.bdt_type

print dy_option
print bdt_option
mcGroups={}
bdts=[]
if dy_option=="amc" : 
  mcGroups = {#'DY': ["DYJetstoLL_HT100","DYJetstoLL_HT100_200","DYJetstoLL_HT200_400","DYJetstoLL_HT400_600","DYJetstoLL_HT600_Inf"],
				'DY_amc': ["DYJetstoLL_amc"],
            'Top': ["TT"],
            'VV': ["WW","WZ","ZZ"],
            'EWKZ': ["EWK_LLJJ"]}
else : 
  mcGroups = {'DY_mdg': ["DYJetstoLL_HT100","DYJetstoLL_HT100_200","DYJetstoLL_HT200_400","DYJetstoLL_HT400_600","DYJetstoLL_HT600_Inf"],
				#'DY': ["DYJetstoLL_amc"],
            'Top': ["TT"],
            'VV': ["WW","WZ","ZZ"],
            'EWKZ': ["EWK_LLJJ"]}
if bdt_option=='atanh' : 
  bdts=['atanhBDT']
else : 
  bdts=['BDT']



expected = {}
expectedStr = {}
JESnorm = {}
for type in ["mu","el"]: 
 combineFile = TFile("shapes/ewkZjj_13TeV_shapes_"+type+"_DY_"+dy_option+"_"+ bdt_option +".root","RECREATE")
 combineFilename = "ewkZjj_13TeV_shapes_"+type+"_DY_"+dy_option+"_"+ bdt_option +".root" 
 for process in ["SingleMuon","SingleElectron"] :
  if process=="SingleMuon" and type!="mu" :
    continue
  if process=="SingleElectron" and type!="el" : 
    continue  
  print "prepare " + type + process
  sourceFile13var = TFile(getTreeLocation() + "ewkZjj_13TeV_"+type+".root")
   # for basePlot in ["BDT","atanhBDT"]:
 
 # for basePlot in ["BDT","atanhBDT"]:
  hBkg = TH1F("hBkg","",300,0,3)
  hSig = TH1F("hSig","",300,0,3)
  binBoundaries = numpy.zeros(NBINS+1,dtype=float)
  binBoundaries[0] = 0
  binBoundaries[NBINS] = 3
  foundHighBinEdge = False
  foundLowBinEdge = True
  B_err2_low = 0.
  B_err2_high = 0.

  for basePlot in bdts:
    sourceFile = sourceFile13var
    basePlot_ = basePlot
    sourceFile.cd()
    for systematic in [""]:
      for name, mcs in mcGroups.iteritems():
        plot = basePlot_ +"_"+ type +"_"+ name+systematic
        DYplot = basePlot_ +"_"+ type +"_DY"+systematic
        thisGroup = merge(sourceFile, plot, mcs)
        if systematic=="" :
          if name=="EWKZ" : 
            hSig = thisGroup
        #    print "signal", thisGroup.Integral()
          if name!="EWKZ" : 
            hBkg.Add(thisGroup)
       #     print name, hBkg.Integral()
    hBkg.Add(merge(sourceFile, basePlot_ + "_" + type + "_interference", ["interference"]))
    print "Bkg", hBkg.GetBinCenter(hBkg.FindLastBinAbove(0))
  #  print hBkg.Integral()
    nBinsFine=300
    for ibin in range(1,nBinsFine):
      if not foundLowBinEdge:
        B_low = hBkg.Integral(1,ibin)
        B_err2_low += pow(hBkg.GetBinError(ibin),2) 
        if (B_low > 0 and sqrt(B_err2_low)/B_low < 0.35):
            binBoundaries[1] = hBkg.GetBinLowEdge(ibin+1)
            foundLowBinEdge = True
      if not foundHighBinEdge:
        B_high = hBkg.Integral(nBinsFine-ibin+1, nBinsFine) 
        B_err2_high += pow(hBkg.GetBinError(nBinsFine-ibin+1),2)
        S_high = hSig.Integral(nBinsFine-ibin+1, nBinsFine) 
      #  if (B_high > 0 and sqrt(B_err2_high)/B_high < 0.35):
        if (B_high > 0 and sqrt(B_err2_high)/B_high < 0.35 ):
      #  if (B_high > 0 and B_high > 1 and sqrt(B_err2_high)/B_high < 0.35):
            binBoundaries[NBINS-1] = hBkg.GetBinLowEdge(nBinsFine-ibin+1)
            foundHighBinEdge = True

    highBinNum = hBkg.FindBin(binBoundaries[NBINS-1])
    lowBinNum = -1
    nBinsLeft = nBinsFine - highBinNum
    nBinsRemain = nBinsLeft % (NBINS-2)
    binBoundaries[0]=0
    binBoundaries[1]=hBkg.GetBinLowEdge(nBinsRemain)
    for i in range(2,NBINS-1):
      binBoundaries[i] = binBoundaries[1] + (i-1)*((binBoundaries[NBINS-1] - binBoundaries[1])/(NBINS-2)) 
    print binBoundaries


  for basePlot in bdts:
    sourceFile = sourceFile13var
    basePlot_ = basePlot
    data = getPlot(sourceFile, "", basePlot_ +"_"+ type +"_"+ "data_obs" )
    rebin_factor = data.GetNbinsX()/NBINS
   # data.Rebin(rebin_factor)
    print data.GetNbinsX()
    print process, data.GetBinCenter(data.FindLastBinAbove(0))
  #  for i in range(0, data.GetNbinsX()+1): 
  #    print data.GetBinContent(i)
    data = data.Rebin(NBINS,"",binBoundaries)
  #  for i in range(0, data.GetNbinsX()+1): 
   #   if data.GetBinContent(i)>10E09 : data.SetBinContent(i,0)
   #   print data.GetBinContent(i)
 #   print process, data.GetBinCenter(data.FindLastBinAbove(0))
    print data.Integral()
    combineFile.cd()
    data.Write()

    for systematic in ["","_CMS_ewkzjj_puWeightUp","_CMS_ewkzjj_puWeightDown","_CMS_ewkzjj_LHE_weights_scaleUp","_CMS_ewkzjj_LHE_weights_scaleDown"]:
      for name, mcs in mcGroups.iteritems():
        if name=="VV" and (systematic=="_CMS_ewkzjj_LHE_weights_scaleUp" or systematic=="_CMS_ewkzjj_LHE_weights_scaleDown") : continue
        plot = basePlot_ +"_"+ type +"_"+ name+systematic
        DYplot = basePlot_ +"_"+ type +"_DY"+systematic
        thisGroup = merge(sourceFile, plot, mcs)
        thisGroup = merge(sourceFile, plot, mcs)
        combineFile.cd()
        rebin_factor = thisGroup.GetNbinsX()/NBINS
       # thisGroup.Rebin(rebin_factor)
        thisGroup = thisGroup.Rebin(NBINS,"",binBoundaries)
#        if systematic=="" : print name, thisGroup.GetBinCenter(thisGroup.FindLastBinAbove(0))
        if systematic=="_CMS_ewkzjj_LHE_weights_scaleUp":
          plot =  basePlot_ +"_"+ type +"_"+ name+"_CMS_ewkzjj_LHE_weights_scale_"+name+"Up"
          DYplot =  basePlot_ +"_"+ type +"_DY"+ "_CMS_ewkzjj_LHE_weights_scale_DY"+"Up"
        if systematic=="_CMS_ewkzjj_LHE_weights_scaleDown":
          plot =  basePlot_ +"_"+ type +"_"+ name+"_CMS_ewkzjj_LHE_weights_scale_"+name+"Down"
          DYplot =  basePlot_ +"_"+ type +"_DY"+ "_CMS_ewkzjj_LHE_weights_scale_DY"+"Down"
        if name.find("DY_")==-1 : thisGroup.Write(plot)
        if name.find("DY_")!=-1 : thisGroup.Write(DYplot)
        if systematic=="" : 
          expectedStr[name] = ('%.3f' % thisGroup.Integral())
          expected[name] = thisGroup.Integral()
        if systematic == "" and name.find("DY")!=-1:
          for bin in range(1, data.GetNbinsX()+1):
         #   print thisGroup.GetBinContent(bin) 
            plotWithBinErrorUp = addErrorToHist(thisGroup,1,bin) 
            rebin_factor = plotWithBinErrorUp.GetNbinsX()/NBINS
         #   plotWithBinErrorUp.Rebin(rebin_factor)
            plotWithBinErrorUp =  plotWithBinErrorUp.Rebin(NBINS,"",binBoundaries)
          #  print plotWithBinErrorUp.GetBinContent(bin) 
            plotWithBinErrorUp.Write(DYplot +"_CMS_ewkzjj_stats_DY_"+type+ "_b" + str(bin) + "Up")
            plotWithBinErrorDown = addErrorToHist(thisGroup,-1,bin)
            rebin_factor = plotWithBinErrorDown.GetNbinsX()/NBINS
         #   plotWithBinErrorDown.Rebin(rebin_factor)
            plotWithBinErrorDown =  plotWithBinErrorDown.Rebin(NBINS,"",binBoundaries)
          #  print plotWithBinErrorDown.GetBinContent(bin) 
            plotWithBinErrorDown.Write(DYplot + "_CMS_ewkzjj_stats_DY_"+type+"_b" + str(bin) + "Down")

 
    interferenceHist = merge(sourceFile, basePlot_ + "_" + type + "_interference", ["interference"])
    rebin_factor = interferenceHist.GetNbinsX()/NBINS
   # interferenceHist.Rebin(rebin_factor)
    interferenceHist = interferenceHist.Rebin(NBINS,"",binBoundaries)
    print NBINS, interferenceHist.GetNbinsX()
    interferenceHist.Write(basePlot_+"_"+type+"_interference")
    expectedStr["interference"] = ('%.3f' % interferenceHist.Integral())
    expected["interference"] = interferenceHist.Integral()
    for intOption in ["","_interference"]:
  #  for intOption in [""]:
      with open("cards/ewkZjj_13TeV_datacard"+ intOption + "_" + type + "_DY"+dy_option+"_"+ bdt_option + ".txt", "wt") as card:
        with open("ewkZjj_all_template" + intOption + ".txt", "rt") as template:
          for line in template:
            line = line.replace('$DIRECTORY', directoryName).replace('$data', '%d'%data.Integral()).replace('$channel',type).replace('$shapesfilename',combineFilename).replace('$bdtname',bdts[0])
            for name, n in expectedStr.iteritems():
              name_used = name
              if name.find('DY')!=-1 : name_used = 'DY'
              line = line.replace(('$'+name_used).ljust(16), n.ljust(16))
            card.write(line)

exit()
