#source $VO_CMS_SW_DIR/cmsset_default.sh
# shopt -s expand_aliases is needed if you want to use the alias 'cmsenv' created by $VO_CMS_SW_DIR/cmsset_default.sh instead of the less mnemonic eval `scramv1 runtime -sh`

source $VO_CMS_SW_DIR/cmsset_default.sh
source /swshare/psit3/etc/profile.d/cms_ui_env.sh  # for bash

export MYCMSENVDIR=/mnt/t3nfs01/data01/shome/nchernya/CMSSW_8_0_19/src
cd $MYCMSENVDIR
eval `scramv1 runtime -sh`
shopt -s expand_aliases 
cmsenv
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib64/dcap 

export MYBATCHDIR=/mnt/t3nfs01/data01/shome/nchernya/VBFZll/main_mva/application/
cd $MYBATCHDIR

dataset_type=(mu el)
root=.root

./run_all2 $1 $TMPDIR/$2 $3


echo $1 $TMPDIR/$2 $3  


#xrdfs t3dcachedb03.psi.ch rm /pnfs/psi.ch/cms/trivcat//store/user/nchernya/VBFZll/mva_v25_reskim/$2_$3$root

xrdcp -f $TMPDIR/$2_$3$root root://t3dcachedb.psi.ch:1094//pnfs/psi.ch/cms/trivcat//store/user/nchernya//VBFZll/mva_v25_reskim/
#xrdcp -f $TMPDIR/$2_$3$root root://t3dcachedb.psi.ch:1094//pnfs/psi.ch/cms/trivcat//store/user/nchernya//VBFZll/skim_for_VBFHmumu/
#rm  $TMPDIR/$2${dataset_type[ $3 ]}$root 


#$ -o /mnt/t3nfs01/data01/shome/nchernya/Hbb/skim_trees/batch_logs2/
#$ -e /mnt/t3nfs01/data01/shome/nchernya/Hbb/skim_trees/batch_logs2/
