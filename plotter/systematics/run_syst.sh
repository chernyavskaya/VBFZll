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
path=dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/mva_syst/
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
["SingleMuon"]=SingleMuon
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
#["interference"]=EWK_LLJJ
["WJetsToLNu"]=WJetsToLNu_madgraph
["ST_tW"]=ST_tW_5f_inclusiveDecays_powheg
["ST_s-channel"]=ST_s-channel_4f_leptonDecays_amc
["ST_t-channel_top_4f_inclusiveDecays"]=ST_t-channel_top_4f_inclusiveDecay_powheg
["ST_t-channel_antitop_4f_inclusiveDecays"]=ST_t-channel_antitop_4f_inclusiveDecays_powheg
)
prefix='main_mva_v24_'
postfix='systematics'
v=v24
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
		#	qsub -q short.q batch.sh  $f ${key} ${region[$current_region]} $data $v $postfix
			echo  $f ${key} ${region[$current_region]} $data $v $postfix
		fi 
		if [ $key != SingleMuon ] && [ $key != SingleMuonB ] && [ $key != SingleMuonC ] && [ $key != SingleMuonD ] && [ $key != SingleMuonE ] && [ $key != SingleMuonF ] && [ $key != SingleMuonG ] && [ $key != SingleElectron ] && [ $key != SingleElectronB ] && [ $key != SingleElectronC ] && [ $key != SingleElectronD ] && [ $key != SingleElectronE ] && [ $key != SingleElectronF ] && [ $key != SingleElectronG ]    
		then
		#	f=$path$prefix${file_names[${key}]}_$v.root
			f=$path$prefix${file_names[${key}]}_${region[$current_region]}.root
	#		qsub -q short.q batch.sh  $f ${key} ${region[$current_region]} $data   $v $postfix
			echo  $f ${key} ${region[$current_region]} $data $v $postfix
		fi
		current_region=$(( $current_region + 1 ))
#		break
	done
#	break
done
