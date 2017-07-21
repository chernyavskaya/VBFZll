export WORKDIR=`pwd`
cd $WORKDIR


declare -A file_names

path=root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat//store/user/nchernya/VBFZll/plotterOutput/v25/
file_name=(
#SingleMuon
#SingleMuonB
#SingleMuonC
#SingleMuonD
#SingleMuonE
#SingleMuonF
#SingleMuonG
#SingleElectron
#SingleElectronB
#SingleElectronC
#SingleElectronD
#SingleElectronE
#SingleElectronF
#SingleElectronG
#EWK_LLJJ
#EWK_LLJJ_herwig
#DYJetstoLL
#DYJetstoLL_amc
#TT
#WW
#WZ
#ZZ
#DYJetstoLL_HT
#DYJetstoLL
DYJetstoLL_amc_0J
DYJetstoLL_amc_1J
DYJetstoLL_amc_2J
#DYJetstoLL_HT100
#DYJetstoLL_HT100_200
#DYJetstoLL_HT200_400
#DYJetstoLL_HT400_600
#DYJetstoLL_HT600_800
#DYJetstoLL_HT800_1200
#DYJetstoLL_HT1200_2500
#DYJetstoLL_HT2500_Inf
#DYJetstoLL_HT600_Inf
#DYJetstoLL_Pt-100_amc
#DYJetstoLL_Pt-100To250_amc
#DYJetstoLL_Pt-250To400_amc
#DYJetstoLL_Pt-400To650_amc
#DYJetstoLL_Pt-650ToInf_amc
ST_tW_top
ST_tW_antitop
ST_s-channel
ST_t-channel_top_4f_inclusiveDecays
ST_t-channel_antitop_4f_inclusiveDecays
#WJetsToLNu_HT100To200
#WJetsToLNu_HT200To400  
#WJetsToLNu_HT400To600
#WJetsToLNu_HT600To800
#WJetsToLNu_HT800To1200
#WJetsToLNu_HT1200To2500
#WJetsToLNu_HT2500ToInf
#WJetsToLNu
#WJetsToLNu_HT
)

prefix=''
postfix=v25_bdt_alldata4_qglweightsnorm_vetoeff_newGapAct_reminiaod
v=v25
ROOT=.root
dataset_type=(mu el)
applyTrigWeight=1
region=(mu el)
output_dir=$TMPDIR
pathhome=/afs/cern.ch/work/n/nchernya/VBFZll/plotter/output_root/


current_region=0
current_file=0
current_trigWeight=0
max_files=8 #16 #35
#	tail -n+3 /afs/cern.ch/work/n/nchernya/VBFZll/plotter/output_txt/${file_name[$current_file]}_${region[$current_region]}_$postfix.txt  | head -n $((1))
#	xrdcp $path${file_name[$current_file]}_${region[$current_region]}_JESnom_$postfix.txt /afs/cern.ch/work/n/nchernya/VBFZll/plotter/output_txt/
while [ $current_file -lt $max_files  ]
do	 
#	xrdcp $path${file_name[$current_file]}_${region[$current_region]}_JESnom_$postfix.txt /afs/cern.ch/work/n/nchernya/VBFZll/plotter/output_txt/
#	xrdcp $path${file_name[$current_file]}_${region[$current_region]}_$postfix.txt /afs/cern.ch/work/n/nchernya/VBFZll/plotter/output_txt/
	echo ${file_name[$current_file]}_${region[$current_region]}
	xrdcp $path${file_name[$current_file]}_${region[$current_region]}_QCDScalenom_JESnom_$postfix.root /afs/cern.ch/work/n/nchernya/VBFZll/plotter/output_root/
	xrdcp $path${file_name[$current_file]}_${region[$current_region]}_QCDScalenom_JESup_$postfix.root /afs/cern.ch/work/n/nchernya/VBFZll/plotter/output_root/
	xrdcp $path${file_name[$current_file]}_${region[$current_region]}_QCDScalenom_JESdown_$postfix.root /afs/cern.ch/work/n/nchernya/VBFZll/plotter/output_root/
	xrdcp $path${file_name[$current_file]}_${region[$current_region]}_QCDScaleup_JESnom_$postfix.root /afs/cern.ch/work/n/nchernya/VBFZll/plotter/output_root/
	xrdcp $path${file_name[$current_file]}_${region[$current_region]}_QCDScaledown_JESnom_$postfix.root /afs/cern.ch/work/n/nchernya/VBFZll/plotter/output_root/
#	tail -n+4 /afs/cern.ch/work/n/nchernya/VBFZll/plotter/output_txt/${file_name[$current_file]}_${region[$current_region]}_QCDScalenom_JESnom_$postfix.txt  | head -n $((1))
#	tail -n+4 /afs/cern.ch/work/n/nchernya/VBFZll/plotter/output_txt/${file_name[$current_file]}_${region[$current_region]}_$postfix.txt  | head -n $((1))
#	xrdcp $path${file_name[$current_file]}_${region[$current_region]}_JESdown_$postfix.root /afs/cern.ch/work/n/nchernya/VBFZll/plotter/output_root/


#	echo
#	echo

	current_file=$((current_file + 1))
done

#WJetsToLNu_HT
list="
DYJetstoLL_amc
ST
"
#DYJetstoLL_HT
#DYJetstoLL_amc

for arg in $list
do

hadd ${pathhome}${arg}_${dataset_type[$current_region]}_QCDScaledown_JESnom_$postfix.root ${pathhome}$arg*_${dataset_type[$current_region]}_QCDScaledown_JESnom_$postfix.root
hadd ${pathhome}${arg}_${dataset_type[$current_region]}_QCDScaleup_JESnom_$postfix.root ${pathhome}$arg*_${dataset_type[$current_region]}_QCDScaleup_JESnom_$postfix.root
hadd ${pathhome}${arg}_${dataset_type[$current_region]}_QCDScalenom_JESnom_$postfix.root ${pathhome}$arg*_${dataset_type[$current_region]}_QCDScalenom_JESnom_$postfix.root
hadd ${pathhome}${arg}_${dataset_type[$current_region]}_QCDScalenom_JESup_$postfix.root ${pathhome}$arg*_${dataset_type[$current_region]}_QCDScalenom_JESup_$postfix.root
hadd ${pathhome}${arg}_${dataset_type[$current_region]}_QCDScalenom_JESdown_$postfix.root ${pathhome}$arg*_${dataset_type[$current_region]}_QCDScalenom_JESdown_$postfix.root

xrdfs t3dcachedb03.psi.ch rm //pnfs/psi.ch/cms/trivcat/store/user//nchernya/VBFZll/plotterOutput/v25/${arg}_${dataset_type[$current_region]}_QCDScaledown_JESnom_$postfix.root
xrdfs t3dcachedb03.psi.ch rm //pnfs/psi.ch/cms/trivcat/store/user//nchernya/VBFZll/plotterOutput/v25/${arg}_${dataset_type[$current_region]}_QCDScaleup_JESnom_$postfix.root
xrdfs t3dcachedb03.psi.ch rm //pnfs/psi.ch/cms/trivcat/store/user//nchernya/VBFZll/plotterOutput/v25/${arg}_${dataset_type[$current_region]}_QCDScalenom_JESnom_$postfix.root
xrdfs t3dcachedb03.psi.ch rm //pnfs/psi.ch/cms/trivcat/store/user//nchernya/VBFZll/plotterOutput/v25/${arg}_${dataset_type[$current_region]}_QCDScalenom_JESup_$postfix.root
xrdfs t3dcachedb03.psi.ch rm //pnfs/psi.ch/cms/trivcat/store/user//nchernya/VBFZll/plotterOutput/v25/${arg}_${dataset_type[$current_region]}_QCDScalenom_JESdown_$postfix.root


#path2=root://t3se01.psi.ch///
path2=root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/

xrdcp -f ${pathhome}${arg}_${dataset_type[$current_region]}_QCDScaledown_JESnom_$postfix.root $path2/store/user/nchernya/VBFZll/plotterOutput/v25/
xrdcp -f ${pathhome}${arg}_${dataset_type[$current_region]}_QCDScaleup_JESnom_$postfix.root $path2/store/user/nchernya/VBFZll/plotterOutput/v25/
xrdcp -f ${pathhome}${arg}_${dataset_type[$current_region]}_QCDScalenom_JESnom_$postfix.root $path2/store/user/nchernya/VBFZll/plotterOutput/v25/
xrdcp -f ${pathhome}${arg}_${dataset_type[$current_region]}_QCDScalenom_JESup_$postfix.root $path2/store/user/nchernya/VBFZll/plotterOutput/v25/
xrdcp -f ${pathhome}${arg}_${dataset_type[$current_region]}_QCDScalenom_JESdown_$postfix.root $path2/store/user/nchernya/VBFZll/plotterOutput/v25/

done
