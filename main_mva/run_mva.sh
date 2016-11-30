export WORKDIR=`pwd`
cd $WORKDIR

g++ run_create_main_tmva_all.C -g -o run_all `root-config --cflags --glibs` 

max_samples_num=2
path=dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat///store/user/nchernya/VBFZll/skimmed/
input_dir=( EWK_LLJJ DYJetstoLL )
ROOT=.root
v24=_v24


current_sample=0
while [ $current_sample -lt $max_samples_num ]
do	
	./run_all $path${input_dir[ $current_sample ]}$v24$ROOT ${input_dir[ $current_sample ]} mu 
	./run_all $path${input_dir[ $current_sample ]}$v24$ROOT ${input_dir[ $current_sample ]} el 

	current_sample=$(( $current_sample + 1 ))
done
