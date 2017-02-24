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
NBINS=20

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
				'DY_amc': ["DYJetstoLL_amc_0J", "DYJetstoLL_amc_1J","DYJetstoLL_amc_2J"],
		#		'DY_amc': ["DYJetstoLL_amc"],
            'Top': ["TT","ST_tW_top","ST_tW_antitop","ST_s-channel","ST_t-channel_top_4f_inclusiveDecays","ST_t-channel_antitop_4f_inclusiveDecays"],
          #  'Top': ["TT"],
         #   'VV': ["WW","WZ","ZZ","WJetsToLNu"],
            'VV': ["WW","WZ","ZZ"],
           # 'VV': ["WW","WZ","ZZ"],
            'EWKZ': ["EWK_LLJJ"]}
else : 
  mcGroups = {
			#'DY_mdg': ["DYJetstoLL_HT100","DYJetstoLL_HT100_200","DYJetstoLL_HT200_400","DYJetstoLL_HT400_600","DYJetstoLL_HT600_Inf"],
				'DY_mdg': ["DYJetstoLL_HT100","DYJetstoLL_HT100_200","DYJetstoLL_HT200_400","DYJetstoLL_HT400_600","DYJetstoLL_HT600_800","DYJetstoLL_HT800_1200","DYJetstoLL_HT1200_2500","DYJetstoLL_HT2500_Inf"],
				#'DY': ["DYJetstoLL_amc"],
            'Top': ["TT","ST_tW_top","ST_tW_antitop","ST_s-channel","ST_t-channel_top_4f_inclusiveDecays","ST_t-channel_antitop_4f_inclusiveDecays"],
            'VV': ["WW","WZ","ZZ"],
           # 'VV': ["WW","WZ","ZZ","WJetsToLNu"],
            'EWKZ': ["EWK_LLJJ"]}
if bdt_option=='atanh' : 
  bdts=['atanhBDT']
else : 
  bdts=['BDT']



expected = {}
expectedStr = {}
JESnorm = {}
for type in ["mu","el"]:
 if type=="mu" : sigAccPS = "1.024"
 if type=="el" : sigAccPS = "1.068"
 combineFile = TFile("shapes/ewkZjj_13TeV_shapes_"+type+"_DY_"+dy_option+"_"+ bdt_option +"_v25alldata.root","RECREATE")
 combineFilename = "ewkZjj_13TeV_shapes_"+type+"_DY_"+dy_option+"_"+ bdt_option +"_v25alldata.root" 
 for process in ["SingleMuon","SingleElectron"] :
  if process=="SingleMuon" and type!="mu" :
    continue
  if process=="SingleElectron" and type!="el" : 
    continue  
  print "prepare " + type + process
  sourceFile13var = TFile(getTreeLocation() + "ewkZjj_13TeV_"+type+"_v25alldata_range053.root")
   # for basePlot in ["BDT","atanhBDT"]:
 
 # for basePlot in ["BDT","atanhBDT"]:
  hBkg = TH1F("hBkg","",750,0.5,3)
  hSig = TH1F("hSig","",750,0.5,3)
  nBinsFine=750
  binBoundaries = numpy.zeros(NBINS+1,dtype=float)
  binBoundaries[0] = 0.5
  binBoundaries[NBINS] = 3
  foundHighBinEdge = False
  foundLowBinEdge = True
  B_err2_low = 0.
  B_err2_high = 0.

  for basePlot in bdts:
    sourceFile = sourceFile13var
    basePlot_ = basePlot
    sourceFile.cd()
    data = getPlot(sourceFile, "", basePlot_ +"_"+ type +"_"+ "data_obs" )
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
            print name, hBkg.Integral()
    hBkg.Add(merge(sourceFile, basePlot_ + "_" + type + "_interference", ["interference"]))
 #   print "Bkg", hBkg.GetBinCenter(hBkg.FindLastBinAbove(0))
    print hBkg.Integral()
    for ibin in range(1,nBinsFine):
      if not foundLowBinEdge:
        B_low = hBkg.Integral(1,ibin)
        B_err2_low += pow(hBkg.GetBinError(ibin),2) 
        if (B_low > 0. and sqrt(B_err2_low)/B_low < 0.35):
            binBoundaries[1] = hBkg.GetBinLowEdge(ibin+1)
            foundLowBinEdge = True
      if not foundHighBinEdge:
        B_high = hBkg.Integral(nBinsFine-ibin+1, nBinsFine) 
        B_err2_high += pow(hBkg.GetBinError(nBinsFine-ibin+1),2)
        S_high = hSig.Integral(nBinsFine-ibin+1, nBinsFine) 
        data_high = data.Integral(nBinsFine-ibin+1, nBinsFine) 
       # if (B_high > 0 and data_high > 1 and sqrt(B_err2_high)/B_high < 0.35):
        if ((S_high+B_high > 0)  and sqrt(B_err2_high)/B_high < 0.35):
       # if ((S_high+B_high) > 1. and sqrt(B_err2_high)/B_high < 0.35 ):
      #  if (B_high > 0 and B_high > 1 and sqrt(B_err2_high)/B_high < 0.35):
            binBoundaries[NBINS-1] = hBkg.GetBinLowEdge(nBinsFine-ibin+1)
            foundHighBinEdge = True

  ####  binBoundaries[NBINS-1] = 2.4
    highBinNum = hBkg.FindBin(binBoundaries[NBINS-1])
    lowBinNum = -1
    nBinsLeft = nBinsFine - highBinNum
    nBinsRemain = nBinsLeft % (NBINS-2)
    binBoundaries[0]=0.5
    print nBinsRemain
    binBoundaries[1]=hBkg.GetBinLowEdge(nBinsRemain)
  #####  binBoundaries[1]=0.3
    for i in range(1,NBINS-1):
      binBoundaries[i] = binBoundaries[0] + i*((binBoundaries[NBINS-1] - binBoundaries[0])/(NBINS-1)) 
    print binBoundaries


  for basePlot in bdts:
    sourceFile = sourceFile13var
    basePlot_ = basePlot
    data = getPlot(sourceFile, "", basePlot_ +"_"+ type +"_"+ "data_obs" )
    rebin_factor = data.GetNbinsX()/NBINS
   # data.Rebin(rebin_factor)
    print data.GetNbinsX()
  #  for i in range(0, data.GetNbinsX()+1): 
  #    print data.GetBinContent(i)
    data = data.Rebin(NBINS,"",binBoundaries)
    print process,"amount of data in the last bin ", data.GetBinContent(data.GetNbinsX())
  #  for i in range(0, data.GetNbinsX()+1): 
   #   if data.GetBinContent(i)>10E09 : data.SetBinContent(i,0)
   #   print data.GetBinContent(i)
 #   print process, data.GetBinCenter(data.FindLastBinAbove(0))
    print data.Integral()
    combineFile.cd()
    data.Write()

    for systematic in ["","_CMS_ewkzjj_puWeightUp","_CMS_ewkzjj_puWeightDown","_CMS_ewkzjj_LHE_weights_scaleUp","_CMS_ewkzjj_LHE_weights_scaleDown","_CMS_ewkzjj_MDG_NLO_corrUp","_CMS_ewkzjj_MDG_NLO_corrDown","_CMS_ewkzjj_JESUp","_CMS_ewkzjj_JESDown","_CMS_ewkzjj_JERUp","_CMS_ewkzjj_JERDown"]:
  #  for systematic in ["","_CMS_ewkzjj_puWeightUp","_CMS_ewkzjj_puWeightDown","_CMS_ewkzjj_LHE_weights_scaleUp","_CMS_ewkzjj_LHE_weights_scaleDown","_CMS_ewkzjj_JESUp","_CMS_ewkzjj_JESDown","_CMS_ewkzjj_JERUp","_CMS_ewkzjj_JERDown"]:
      for name, mcs in mcGroups.iteritems():
        if (name=="VV" or name=="Top" ) and (systematic=="_CMS_ewkzjj_LHE_weights_scaleUp" or systematic=="_CMS_ewkzjj_LHE_weights_scaleDown") : continue
        if name!="DY_mdg" and (systematic=="_CMS_ewkzjj_MDG_NLO_corrUp" or systematic=="_CMS_ewkzjj_MDG_NLO_corrDown" ) : continue
        if name!="interference" and (systematic=="_CMS_ewkzjj_int_shapeUp" or systematic=="_CMS_ewkzjj_int_shapeDown" ) : continue
        plot = basePlot_ +"_"+ type +"_"+ name+systematic
        DYplot = basePlot_ +"_"+ type +"_DY"+systematic
        thisGroup = merge(sourceFile, plot, mcs)
        combineFile.cd()
       # rebin_factor = thisGroup.GetNbinsX()/NBINS
       # thisGroup.Rebin(rebin_factor)
        thisGroup = thisGroup.Rebin(NBINS,"",binBoundaries)
        if systematic =="" : print "Last bin in process has that many events", name, thisGroup.GetBinContent(thisGroup.GetNbinsX())
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
          print 'number of bins for DY', thisGroup.GetNbinsX()
          for bin in range(1, thisGroup.GetNbinsX()+1):
         #   print thisGroup.GetBinContent(bin) 
            plotWithBinErrorUp = addErrorToHist(thisGroup,1,bin) 
          #  rebin_factor = plotWithBinErrorUp.GetNbinsX()/NBINS
         #   plotWithBinErrorUp.Rebin(rebin_factor)
       #####     plotWithBinErrorUp =  plotWithBinErrorUp.Rebin(NBINS,"",binBoundaries)
          #  print plotWithBinErrorUp.GetBinContent(bin) 
            plotWithBinErrorUp.Write(DYplot +"_CMS_ewkzjj_stats_DY_"+type+ "_b" + str(bin) + "Up")
            plotWithBinErrorDown = addErrorToHist(thisGroup,-1,bin)
         #   rebin_factor = plotWithBinErrorDown.GetNbinsX()/NBINS
         #   plotWithBinErrorDown.Rebin(rebin_factor)
        #####    plotWithBinErrorDown =  plotWithBinErrorDown.Rebin(NBINS,"",binBoundaries)
          #  print plotWithBinErrorDown.GetBinContent(bin) 
            plotWithBinErrorDown.Write(DYplot + "_CMS_ewkzjj_stats_DY_"+type+"_b" + str(bin) + "Down")

 
    interferenceHist = merge(sourceFile, basePlot_ + "_" + type + "_interference", ["interference"])
   # rebin_factor = interferenceHist.GetNbinsX()/NBINS
   # interferenceHist.Rebin(rebin_factor)
    interferenceHist = interferenceHist.Rebin(NBINS,"",binBoundaries)
    print "Last bin in process has that many events in interference", interferenceHist.GetBinContent(interferenceHist.GetNbinsX())
    interferenceHist.Write(basePlot_+"_"+type+"_interference")
    interferenceHistUp = merge(sourceFile, basePlot_ + "_" + type + "_interference_CMS_ewkzjj_int_shapeUp", ["interference"])
    interferenceHistUp = interferenceHistUp.Rebin(NBINS,"",binBoundaries)
    interferenceHistUp.Write(basePlot_+"_"+type+"_interference_CMS_ewkzjj_int_shapeUp")
    interferenceHistDown = merge(sourceFile, basePlot_ + "_" + type + "_interference_CMS_ewkzjj_int_shapeDown", ["interference"])
    interferenceHistDown = interferenceHistDown.Rebin(NBINS,"",binBoundaries)
    interferenceHistDown.Write(basePlot_+"_"+type+"_interference_CMS_ewkzjj_int_shapeDown")
    expectedStr["interference"] = ('%.3f' % interferenceHist.Integral())
    expected["interference"] = interferenceHist.Integral()
    for intOption in ["","_interference"]:
  #  for intOption in [""]:
      with open("cards/ewkZjj_13TeV_datacard"+ intOption + "_" + type + "_DY"+dy_option+"_"+ bdt_option + "_v25alldata.txt", "wt") as card:
        with open("ewkZjj_all_template" + intOption + ".txt", "rt") as template:
          for line in template:
            line = line.replace('$DIRECTORY', directoryName).replace('$data', '%d'%data.Integral()).replace('$channel',type).replace('$shapesfilename',combineFilename).replace('$bdtname',bdts[0]).replace('$SigPS',sigAccPS)
            if line.find("MDG_NLO_corr")!=-1 and dy_option=="amc" : continue
            for name, n in expectedStr.iteritems():
              name_used = name
              if name.find('DY')!=-1 : name_used = 'DY'
              line = line.replace(('$'+name_used).ljust(16), n.ljust(16))
            card.write(line)

exit()
