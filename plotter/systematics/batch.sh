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

export MYBATCHDIR=//mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/systematics/
cd $MYBATCHDIR

root=.root

./syst $1 $2 $3 $4 $5 $6  $TMPDIR
xrdfs t3se01.psi.ch mkdir -p /store/user/nchernya/VBFZll/systematics/
xrdfs t3se01.psi.ch mkdir -p /store/user/nchernya/VBFZll/systematics/$5/

xrdfs t3se01.psi.ch rm /store/user/nchernya//VBFZll/systematics/$5//$2_$3_$5_${6}.root


xrdcp $TMPDIR/$2_$3_$5_${6}.root root://t3se01.psi.ch//store/user/nchernya/VBFZll/systematics/$5/

#$ -o /mnt/t3nfs01/data01/shome/nchernya/Hbb/skim_trees/batch_logs/
#$ -e /mnt/t3nfs01/data01/shome/nchernya/Hbb/skim_trees/batch_logs/
