export WORKDIR=`pwd`
cd $WORKDIR

g++ system.C -g -o syst `root-config --cflags --glibs`  -lMLP -lXMLIO -lTMVA  

declare -A file_names

cp syst /mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/systematics/
cp  batch.sh /mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/systematics/
cp  ../EWcorr.C	/mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/systematics/
#cp ../muon_corrections/rochcor2016.cc 	/mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/systematics/
#cp ../muon_corrections/rochcor2016.h 	/mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/systematics/
#cp ../muon_corrections/RoccoR.cc /mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/systematics/
#cp ../muon_corrections/RoccoR.h /mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/systematics/


#path=dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/skimmed/
path=dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/mva_v25_reskim/
file_names=(
["interference"]=EWK_LL_JJ
["DYJetstoLL_amc_0J"]=DYJetstoLL_amc_0J
["DYJetstoLL_amc_1J"]=DYJetstoLL_amc_1J
["DYJetstoLL_amc_2J"]=DYJetstoLL_amc_2J
#["DYJetstoLL"]=DYJetstoLL_madgraph
#["DYJetstoLL_HT100"]=DYJetstoLL_madgraph
#["DYJetstoLL_HT100_200"]=DYJetstoLL_HT100to200
#["DYJetstoLL_HT200_400"]=DYJetstoLL_HT200to400
#["DYJetstoLL_HT400_600"]=DYJetstoLL_HT400to600
#["DYJetstoLL_HT600_800"]=DYJetstoLL_HT600to800
#["DYJetstoLL_HT800_1200"]=DYJetstoLL_HT800to1200
#["DYJetstoLL_HT1200_2500"]=DYJetstoLL_HT1200to2500
#["DYJetstoLL_HT2500_Inf"]=DYJetstoLL_HT2500toInf
["EWK_LLJJ"]=EWK_LL_JJ
["EWK_LLJJ_herwig"]=EWK_LL_JJ_herwig
["SingleMuon"]=SingleMuon_reminiaod
["SingleElectron"]=SingleElectron_reminiaod
["ST_tW_top"]=ST_tW_top
["ST_tW_antitop"]=ST_tW_antitop
["ST_s-channel"]=ST_s
["ST_t-channel_top_4f_inclusiveDecays"]=ST_t_top
["ST_t-channel_antitop_4f_inclusiveDecays"]=ST_t_antitop
["TT"]=TT
["WW"]=WW
["WZ"]=WZ
["ZZ"]=ZZ
["WJetsToLNu"]=WJetsToLnu_madgraph
#["TTZToLLNuNu"]=TTZToLLNuNu
#["tZq_ll"]=tZq_ll
#["EWKinterference"]=LLJJ_EWK_5f_LO_13TeV



)
prefix='main_mva_v25_'
postfix='systematics_alldata5'
v=v25
ROOT=.root
region=(mu el)
output_dir=$TMPDIR


for key in ${!file_names[@]}; do
	current_region=0
#	while [ $current_region -lt 1  ] 
	while [ $current_region -lt  2 ] 
	do
		data=0
		if [ $key == SingleMuon ] || [ $key == SingleMuonB ] || [ $key == SingleMuonC ] || [ $key == SingleMuonD ] || [ $key == SingleMuonE ] || [ $key == SingleMuonF ] || [ $key == SingleMuonG ] ||[ $key == SingleElectron ] || [ $key == SingleElectron ] || [ $key == SingleElectronB ] || [ $key == SingleElectronC ] || [ $key == SingleElectronD ] || [ $key == SingleElectronE ] || [ $key == SingleElectronF ] || [ $key == SingleElectronG ] 
		then 
	 		data=1
		fi
		if [ $key == SingleMuon ] || [ $key == SingleMuonB ] || [ $key == SingleMuonC ] || [ $key == SingleMuonD ] || [ $key == SingleMuonE ] || [ $key == SingleMuonF ] || [ $key == SingleMuonG ] || [ $key == SingleElectron ]  || [ $key == SingleElectronB ] || [ $key == SingleElectronC ] || [ $key == SingleElectronD ] || [ $key == SingleElectronE ] || [ $key == SingleElectronF ] || [ $key == SingleElectronG ] 
		then
		#	f=$path$prefix${file_names[${key}]}_$v.root
			f=$path$prefix${file_names[${key}]}_${region[$current_region]}.root
			qsub -q short.q batch.sh  $f ${key} ${region[$current_region]} $data $v $postfix
	#		echo  $f ${key} ${region[$current_region]} $data $v $postfix
		fi 
		if [ $key != SingleMuon ] && [ $key != SingleMuonB ] && [ $key != SingleMuonC ] && [ $key != SingleMuonD ] && [ $key != SingleMuonE ] && [ $key != SingleMuonF ] && [ $key != SingleMuonG ] && [ $key != SingleElectron ] && [ $key != SingleElectronB ] && [ $key != SingleElectronC ] && [ $key != SingleElectronD ] && [ $key != SingleElectronE ] && [ $key != SingleElectronF ] && [ $key != SingleElectronG ]    
		then
		#	f=$path$prefix${file_names[${key}]}_$v.root
			f=$path$prefix${file_names[${key}]}_${region[$current_region]}.root
			qsub -q short.q batch.sh  $f ${key} ${region[$current_region]} $data   $v $postfix
	#		echo  $f ${key} ${region[$current_region]} $data $v $postfix
		fi
		current_region=$(( $current_region + 1 ))
#		break
	done
#	break
done
