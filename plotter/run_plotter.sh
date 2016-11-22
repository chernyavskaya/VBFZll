export WORKDIR=`pwd`
cd $WORKDIR

g++ plotter_vbfzll.C -g -o plot `root-config --cflags --glibs`  -lMLP -lXMLIO -lTMVA  

declare -A file_names

cp plot /mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/
cp  batch.sh /mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/
cp  EWcorr.C	/mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/
cp muon_corrections/rochcor2016.cc 	/mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/
cp muon_corrections/rochcor2016.h 	/mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/
cp muon_corrections/RoccoR.cc /mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/
cp muon_corrections/RoccoR.h /mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/


path=dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/nchernya/VBFZll/skimmed/
file_names=(
#["DYJetstoLL"]=DYJetstoLL
["DYJetstoLL_amc"]=DYJetsToLL_amc_full
["DYJetstoLL_HT100"]=DYJetstoLL
["DYJetstoLL_HT100_200"]=DYJetstoLL_HT100_200_full
["DYJetstoLL_HT200_400"]=DYJetstoLL_HT200_400_full
["DYJetstoLL_HT400_600"]=DYJetstoLL_HT400_600_full
["DYJetstoLL_HT600_Inf"]=DYJetstoLL_HT600_Inf_full
#["DYJetstoLL_Pt-100_amc"]=DYJetstoLL_amc
#["DYJetstoLL_Pt-100To250_amc"]=DYJetsToLL_Pt-100To250_amc
#["DYJetstoLL_Pt-250To400_amc"]=DYJetsToLL_Pt-250To400_amc
#["DYJetstoLL_Pt-400To650_amc"]=DYJetsToLL_Pt-400To650_amc
#["DYJetstoLL_Pt-650ToInf_amc"]=DYJetsToLL_Pt-650ToInf_amc
#["TT"]=TT
#["WW"]=WW
#["WZ"]=WZ
#["ZZ"]=ZZ
#["SingleMuon"]=SingleMuon_full
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
)
prefix=''
postfix='ewk_mucorr_PtsumEtaQQ'
v=v24
ROOT=.root
region=(mu el)
applyJESWeight=1
applyQCDWeight=1
JESWeightNom=(nom up down)
QCDWeightNom=(nom up down)
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
		if [ $key == SingleMuon ] || [ $key == SingleMuonB ] || [ $key == SingleMuonC ] || [ $key == SingleMuonD ] || [ $key == SingleMuonE ] || [ $key == SingleMuonF ] || [ $key == SingleMuonG ] || [ $key == SingleElectron ]  || [ $key == SingleElectronB ] || [ $key == SingleElectronC ] || [ $key == SingleElectronD ] || [ $key == SingleElectronE ] || [ $key == SingleElectronF ] || [ $key == SingleElectronG ] || [ $applyJESWeight -eq 0 ]
		then
			f=$path$prefix${file_names[${key}]}_$v.root
			qsub -q short.q batch.sh  $f ${key} ${region[$current_region]} $data 0 nom 0 nom $v $postfix
		#	echo  $f ${key} ${region[$current_region]} $data 0 nom 0 nom $v $postfix
		fi 
		if [ $key != SingleMuon ] && [ $key != SingleMuonB ] && [ $key != SingleMuonC ] && [ $key != SingleMuonD ] && [ $key != SingleMuonE ] && [ $key != SingleMuonF ] && [ $key != SingleMuonG ] && [ $key != SingleElectron ] && [ $key != SingleElectronB ] && [ $key != SingleElectronC ] && [ $key != SingleElectronD ] && [ $key != SingleElectronE ] && [ $key != SingleElectronF ] && [ $key != SingleElectronG ]     && [ $applyJESWeight -eq 1 ]
		then
			current_JESWeight=0
			current_QCDWeight=0
	#		while [ $current_QCDWeight -lt  1 ] 
			while [ $current_QCDWeight -lt  1 ] 
			do
				f=$path$prefix${file_names[${key}]}_$v.root
				qsub -q short.q batch.sh  $f ${key} ${region[$current_region]} $data  $applyQCDWeight ${QCDWeightNom[$current_QCDWeight]} 0 nom $v $postfix
	#			echo  $f ${key} ${region[$current_region]} $data  $applyQCDWeight ${QCDWeightNom[$current_QCDWeight]} 0 nom $v $postfix
				current_QCDWeight=$(( $current_QCDWeight + 1 ))
			done
	#		while [ $current_JESWeight -lt  1 ] 
			while [ $current_JESWeight -lt  1 ] 
			do
				f=$path$prefix${file_names[${key}]}_$v.root
				qsub -q short.q batch.sh  $f ${key} ${region[$current_region]} $data 0 nom $applyJESWeight ${JESWeightNom[$current_JESWeight]} $v $postfix
			#	echo $f ${key} ${region[$current_region]} $data  $applyJESWeight ${JESWeightNom[$current_JESWeight]} $v $postfix
				current_JESWeight=$(( $current_JESWeight + 1 ))
			done
		fi
		current_region=$(( $current_region + 1 ))
#		break
	done
#	break
done
