export WORKDIR=`pwd`
cd $WORKDIR

g++ plotter_vbfzll_pdf.C -g -o pdf `root-config --cflags --glibs`  -lMLP -lXMLIO -lTMVA  

declare -A file_names

cp pdf /mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/
cp  batch.sh /mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/
cp  EWcorr.C	/mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/
#cp muon_corrections/rochcor2016.cc 	/mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/
#cp muon_corrections/rochcor2016.h 	/mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/
#cp muon_corrections/RoccoR.cc /mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/
#cp muon_corrections/RoccoR.h /mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/

#path=dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/mva_v25_axis2jet2q/
path=dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/mva_v25_reskim/
#path=dcap://t3dcachedb.psi.ch:1094//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/mva_v25_new/
file_names=(
#["DYJetstoLL"]=DYJetstoLL
#["DYJetstoLL_amc"]=DYJetsToLL_amc_full
#["DYJetstoLL_HT100"]=DYJetstoLL
#["DYJetstoLL_HT100_200"]=DYJetstoLL_HT100_200_full
#["DYJetstoLL_HT200_400"]=DYJetstoLL_HT200_400_full
#["DYJetstoLL_HT400_600"]=DYJetstoLL_HT400_600_full
#["DYJetstoLL_HT600_Inf"]=DYJetstoLL_HT600_Inf_full
#["DYJetstoLL_Pt-100_amc"]=DYJetstoLL_amc
#["DYJetstoLL_Pt-100To250_amc"]=DYJetsToLL_Pt-100To250_amc
#["DYJetstoLL_Pt-250To400_amc"]=DYJetsToLL_Pt-250To400_amc
#["DYJetstoLL_Pt-400To650_amc"]=DYJetsToLL_Pt-400To650_amc
#["DYJetstoLL_Pt-650ToInf_amc"]=DYJetsToLL_Pt-650ToInf_amc
#["TT"]=TT
#["WW"]=WW
#["WZ"]=WZ
#["ZZ"]=ZZ
#["SingleMuon"]=SingleMuon
#["SingleMuonB"]=SingleMuonB
#["SingleMuonC"]=SingleMuonC
#["SingleMuonD"]=SingleMuonD
#["SingleMuonE"]=SingleMuonE
#["SingleMuonF"]=SingleMuonF
#["SingleMuonG"]=SingleMuonG
#["SingleElectron"]=SingleElectron
#["SingleElectronB"]=SingleElectronB
#["SingleElectronC"]=SingleElectronC
#["SingleElectronD"]=SingleElectronD
#["SingleElectronE"]=SingleElectronE
#["SingleElectronF"]=SingleElectronF
#["SingleElectronG"]=SingleElectronG
#["EWK_LLJJ"]=EWK_LLJJ
#######["interference"]=EWK_LLJJ
###["QCD_HT100to200"]=QCD_HT100to200_full
###["QCD_HT200to300"]=QCD_HT200to300_full
###["QCD_HT300to500"]=QCD_HT300to500_full
###["QCD_HT500to700"]=QCD_HT500to700_full
###["QCD_HT700to1000"]=QCD_HT700to1000_full
###["QCD_HT1000to1500"]=QCD_HT1000to1500_full
###["QCD_HT1500to2000"]=QCD_HT1500to2000_full
###["QCD_HT2000toInf"]=QCD_HT2000toInf_full

#["WJetsToLNu_HT100"]=WJetsToLNu_madgraph 
#["WJetsToLNu_HT100To200"]=WJetsToLNu_HT-100To200_madgraph_full 
#["WJetsToLNu_HT200To400"]=WJetsToLNu_HT-200To400_madgraph_full  
#["WJetsToLNu_HT400To600"]=WJetsToLNu_HT-400To600_madgraph_full
#["WJetsToLNu_HT600To800"]=WJetsToLNu_HT-600To800_madgraph_full
#["WJetsToLNu_HT800To1200"]=WJetsToLNu_HT-800To1200_madgraph_full
#["WJetsToLNu_HT1200To2500"]=WJetsToLNu_HT-1200To2500_madgraph_full
#["WJetsToLNu_HT2500ToInf"]=WJetsToLNu_HT-2500ToInf_madgraph_full
#["WJetsToLNu"]=WJetsToLNu_madgraph
#["WJetsToLNu_amc"]=WJetsToLNu_amcatnlo

#["ST_tW"]=ST_tW_5f_inclusiveDecays_powheg
#["ST_s-channel"]=ST_s-channel_4f_leptonDecays_amc
#["ST_t-channel_top_4f_inclusiveDecays"]=ST_t-channel_top_4f_inclusiveDecay_powheg
#["ST_t-channel_antitop_4f_inclusiveDecays"]=ST_t-channel_antitop_4f_inclusiveDecays_powheg
###["ST_t-channel_top_4f_leptonDecays"]=ST_t-channel_top_4f_leptonDecays_powheg
###["ST_t-channel_antitop_4f_leptonDecays"]=ST_t-channel_antitop_4f_leptonDecays_powheg



#["DYJetstoLL_amc_0J"]=DYJetstoLL_amc_0J
#["DYJetstoLL_amc_1J"]=DYJetstoLL_amc_1J
["DYJetstoLL_amc_2J"]=DYJetstoLL_amc_2J
#["DYJetstoLL_amc_2J"]=DYJetstoLL_amc_2J_all
#["DYJetstoLL"]=DYJetstoLL_madgraph
#["DYJetstoLL_HT100"]=DYJetstoLL_madgraph
#["DYJetstoLL_HT100_200"]=DYJetstoLL_HT100to200
#["DYJetstoLL_HT200_400"]=DYJetstoLL_HT200to400
#["DYJetstoLL_HT400_600"]=DYJetstoLL_HT400to600
#["DYJetstoLL_HT600_800"]=DYJetstoLL_HT600to800
#["DYJetstoLL_HT800_1200"]=DYJetstoLL_HT800to1200
#["DYJetstoLL_HT1200_2500"]=DYJetstoLL_HT1200to2500
#["DYJetstoLL_HT2500_Inf"]=DYJetstoLL_HT2500toInf
#["EWK_LLJJ"]=EWK_LL_JJ
#["EWK_LLJJ_herwig"]=EWK_LL_JJ_herwig
#["SingleMuon"]=SingleMuon_reminiaod
#["SingleElectron"]=SingleElectron_reminiaod
#["ST_tW_top"]=ST_tW_top
#["ST_tW_antitop"]=ST_tW_antitop
#["ST_s-channel"]=ST_s
#["ST_t-channel_top_4f_inclusiveDecays"]=ST_t_top
#["ST_t-channel_antitop_4f_inclusiveDecays"]=ST_t_antitop
#["TT"]=TT
#["WW"]=WW
#["WZ"]=WZ
#["ZZ"]=ZZ
#["WJetsToLNu"]=WJetsToLnu_madgraph
#["interference"]=EWK_LL_JJ
#["TTZToLLNuNu"]=TTZToLLNuNu
#["tZq_ll"]=tZq_ll


)
prefix='main_mva_v25_'
#postfix='ewk_mucorr_MqqLog_bdt'
#postfix='ewk_mucorr_nocorr_bdt_oldxsec'
postfix='_pdf_unc2try_reminiaod'
v='v25'
ROOT=.root
region=(mu el)
applyJESWeight=0
applyQCDWeight=0
JESWeightNom=(nom up down)
QCDWeightNom=(nom up down)
output_dir=$TMPDIR


for key in ${!file_names[@]}; do
	current_region=0
#	while [ $current_region -lt 1  ] 
	while [ $current_region -lt  1 ] 
	do
		data=0
		if [ $key == SingleMuon ] || [ $key == SingleMuonB ] || [ $key == SingleMuonC ] || [ $key == SingleMuonD ] || [ $key == SingleMuonE ] || [ $key == SingleMuonF ] || [ $key == SingleMuonG ] ||[ $key == SingleElectron ] || [ $key == SingleElectron ] || [ $key == SingleElectronB ] || [ $key == SingleElectronC ] || [ $key == SingleElectronD ] || [ $key == SingleElectronE ] || [ $key == SingleElectronF ] || [ $key == SingleElectronG ] 
		then 
	 		data=1
		fi
		if [ $key == SingleMuon ] || [ $key == SingleMuonB ] || [ $key == SingleMuonC ] || [ $key == SingleMuonD ] || [ $key == SingleMuonE ] || [ $key == SingleMuonF ] || [ $key == SingleMuonG ] || [ $key == SingleElectron ]  || [ $key == SingleElectronB ] || [ $key == SingleElectronC ] || [ $key == SingleElectronD ] || [ $key == SingleElectronE ] || [ $key == SingleElectronF ] || [ $key == SingleElectronG ] || [ $applyJESWeight -eq 0 ]
		then
		#	f=$path$prefix${file_names[${key}]}_$v.root
			f=$path$prefix${file_names[${key}]}_${region[$current_region]}.root
			qsub -q all.q batch_pdf.sh  $f ${key} ${region[$current_region]} $data  $v $postfix
			echo  $f ${key} ${region[$current_region]} $data  $v $postfix
		fi 
		current_region=$(( $current_region + 1 ))
#		break
	done
#	break
done
