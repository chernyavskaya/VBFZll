export WORKDIR=`pwd`
cd $WORKDIR

g++ run_make_all2.C -g -o run_all2 `root-config --cflags --glibs`  -lMLP -lXMLIO -lTMVA  

cp run_all2 /mnt/t3nfs01/data01/shome/nchernya/VBFZll/main_mva/application/

input_dir=(
#WJetsToLnu_madgraph
#SingleMuon_reminiaod
#SingleElectron_reminiaod
#EWK_LL_JJ
#EWK_LL_JJ_herwig
#DYJetstoLL_madgraph
#DYJetstoLL_HT100to200
#DYJetstoLL_HT200to400
#DYJetstoLL_HT400to600
#DYJetstoLL_HT600to800
#DYJetstoLL_HT800to1200
#DYJetstoLL_HT1200to2500
#DYJetstoLL_HT2500toInf
#DYJetstoLL_amc_0J
#DYJetstoLL_amc_1J
#DYJetstoLL_amc_2J_ext
#DYJetstoLL_amc_2J
#ZZ
#WZ
#WW
#TT
#ST_t_top
#ST_tW_antitop
#ST_s
#ST_tW_top
#ST_t_antitop
###TTZToLLNuNu
####tZq_ll
#GluGlu_HToMuMu_passall
#GluGlu_HToMuMu
#VBF_HToMuMu_passall
#VBF_HToMuMu
#LLJJ_EWK_5f_LO_13TeV
EWK_LLJJ_aTGC_full2
)
max_samples_num=1 # 26   #25 #26  #42  #42  #42  #14

#path=dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat///store//user/nchernya/VBFZll/skimmed/
path=dcap://t3se01.psi.ch:22125///pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/v25_reskim/
#path=dcap://t3se01.psi.ch:22125///pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/v25/
ROOT=.root
file_end=_v25_reskim
#file_end=_v25
slash=/
mva=main_mva_v25_

current_sample=0
while [ $current_sample -lt $max_samples_num ]
#while [ $current_sample -lt 2 ]
do	
#	echo $path${input_dir[ $current_sample ]}$file_end$ROOT $mva${input_dir[ $current_sample ]} mu
	qsub -q short.q batch2_2.sh $path${input_dir[ $current_sample ]}$file_end$ROOT $mva${input_dir[ $current_sample ]} mu
	qsub -q short.q batch2_2.sh $path${input_dir[ $current_sample ]}$file_end$ROOT $mva${input_dir[ $current_sample ]} el
	current_sample=$(( $current_sample + 1 ))
done

