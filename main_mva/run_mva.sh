export WORKDIR=`pwd`
cd $WORKDIR

g++ run_create_main_tmva_all.C -g -o run_all `root-config --cflags --glibs` 

max_samples_num=1 #9
path=dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat///store/user/nchernya/VBFZll/v25/
input_dir=( EWK_LL_JJ DYJetstoLL_madgraph DYJetstoLL_HT100to200 DYJetstoLL_HT200to400 DYJetstoLL_HT400to600 DYJetstoLL_HT600to800 DYJetstoLL_HT800to1200 DYJetstoLL_HT1200to2500 DYJetstoLL_HT2500toInf )
ROOT=.root
v=_v25


current_sample=0
while [ $current_sample -lt $max_samples_num ]
do
	file=$path${input_dir[ $current_sample ]}$v$ROOT
#	./run_all $file ${input_dir[ $current_sample ]} mu 
#	./run_all $file  ${input_dir[ $current_sample ]} el test
	./run_all $file  ${input_dir[ $current_sample ]} el train
#	./run_all $file  ${input_dir[ $current_sample ]} mu test
#	./run_all $file  ${input_dir[ $current_sample ]} mu train

	current_sample=$(( $current_sample + 1 ))
done
