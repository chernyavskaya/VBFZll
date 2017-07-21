import ROOT
from ROOT import gROOT
from ROOT import gStyle
import sys


channel = sys.argv[1]
vtypeSim = int(sys.argv[2])

f0 = ROOT.TFile.Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat//store/user/nchernya/VBFZll/mva_v25_reskim/main_mva_v25_DYJetstoLL_amc_0J_%s.root"%channel)
f1 = ROOT.TFile.Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat//store/user/nchernya/VBFZll/mva_v25_reskim/main_mva_v25_DYJetstoLL_amc_1J_%s.root"%channel)
f2 = ROOT.TFile.Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat//store/user/nchernya/VBFZll/mva_v25_reskim/main_mva_v25_DYJetstoLL_amc_2J_%s.root"%channel)

list_nom=[]
list_denom=[]
list_integrals=[0.278555,1.31182,9.89926]
if channel=='el' :
	list_integrals=[0.0977376,0.486609,3.85313]
	
for i in range(0,len(list_integrals)):
	list_integrals[i]=list_integrals[i]*35900
#	print list_integrals[i]


list_nom.append((f0.Get("tree")).GetEntries("VtypeSim==%i&&PassSelection_nom==1"%vtypeSim)*1.0)
list_denom.append((f0.Get("tree")).GetEntries("PassSelection_nom==1")*1.0)
list_nom.append((f1.Get("tree")).GetEntries("VtypeSim==%i&&PassSelection_nom==1"%vtypeSim)*1.0)
list_denom.append((f1.Get("tree")).GetEntries("PassSelection_nom==1")*1.0)
list_nom.append((f2.Get("tree")).GetEntries("VtypeSim==%i&&PassSelection_nom==1"%vtypeSim)*1.0)
list_denom.append((f2.Get("tree")).GetEntries("PassSelection_nom==1")*1.0)

sum_tau=0.
sum=0.
for i in range(0, len(list_nom)):
#	print list_nom[i],list_denom[i],list_integrals[i]
	sum_tau=sum_tau + (list_nom[i]/list_denom[i]*list_integrals[i])
for i in range(0, len(list_nom)):
	sum=sum+list_integrals[i]

print sum_tau, sum_tau/sum

