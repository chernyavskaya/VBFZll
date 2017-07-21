export WORKDIR=`pwd`
cd $WORKDIR
g++ tmva.C -g -o run_Nm1_tmva  -lTMVA `root-config --cflags --glibs`
cp run_Nm1_tmva /mnt/t3nfs01/data01/shome/nchernya/VBFZll/main_mva  
cp TMVAClassification_main.C /mnt/t3nfs01/data01/shome/nchernya/VBFZll/main_mva
cp tmva.C /mnt/t3nfs01/data01/shome/nchernya/VBFZll/main_mva
cp TMVAGui.C /mnt/t3nfs01/data01/shome/nchernya/VBFZll/main_mva 
cp tmvaglob.C  /mnt/t3nfs01/data01/shome/nchernya/VBFZll/main_mva
cp  efficiencies.C  /mnt/t3nfs01/data01/shome/nchernya/VBFZll/main_mva 


max=1 #7 #8
#variables_names=( Mqq DeltaEtaQQ Jet2q_pt Jet1q_pt Jet1q_leadTrackPt Jet2q_leadTrackPt axis2_jet1 axis2_jet2 qq_pt Jet3_pt RptHard Zll_pt Zll_zstar all ) #1
#variables_names=( Mqq DeltaEtaQQ Jet2q_pt Jet1q_pt Jet1q_leadTrackPt Jet2q_leadTrackPt axis2_jet1 axis2_jet2 qq_pt RptHard Zll_pt Zll_zstar all ) #1
#variables_names=( Mqq DeltaEtaQQ Jet2q_pt axis2_jet1 axis2_jet2 qq_pt RptHard Zll_zstar all ) #1
#variables_names=( Mqq DeltaEtaQQ axis2_jet1 axis2_jet2 qq_pt RptHard Zll_zstar all ) #1
#variables_names=( Mqq DeltaEtaQQ axis2_jet1 axis2_jet2 qq_pt Zll_zstar all ) #1
#variables_names=(all Mqq DeltaEtaQQ Jet2q_pt axis2_jet1 axis2_jet2 qq_pt Zll_zstar ) #1
variables_names=(Zll_zstar) #1
counter=0
while  [ $counter -lt $max ]
do
	nTrees=200
	MinNodeSize=7 # it is in % already
	maxDepth=3
	echo $nTrees  $MinNodeSize $maxDepth
	qsub -q short.q batch2.sh ${variables_names[counter]} elmu $nTrees $MinNodeSize $maxDepth 
	nTrees=200
	MinNodeSize=20   # it is in % already
	maxDepth=3
	echo $nTrees  $MinNodeSize $maxDepth
#	qsub -q short.q batch2.sh ${variables_names[counter]} el $nTrees $MinNodeSize $maxDepth
	counter=$(( $counter + 1 ))	
done
