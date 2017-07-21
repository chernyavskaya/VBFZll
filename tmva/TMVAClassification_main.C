// @(#)root/tmva $Id$
/**********************************************************************************
 * Project   : TMVA - a ROOT-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Root Macro: TMVAClassification                                                 *
 *                                                                                *
 * This macro provides examples for the training and testing of the               *
 * TMVA classifiers.                                                              *
 *                                                                                *
 * As input data is used a toy-MC sample consisting of four Gaussian-distributed  *
 * and linearly correlated input variables.                                       *
 *                                                                                *
 * The methods to be used can be switched on and off by means of booleans, or     *
 * via the prompt command, for example:                                           *
 *                                                                                *
 *    root -l ./TMVAClassification.C\(\"Fisher,Likelihood\"\)                     *
 *                                                                                *
 * (note that the backslashes are mandatory)                                      *
 * If no method given, a default set of classifiers is used.                      *
 *                                                                                *
 * The output file "TMVA.root" can be analysed with the use of dedicated          *
 * macros (simply say: root -l <macro.C>), which can be conveniently              *
 * invoked through a GUI that will appear at the end of the run of this macro.    *
 * Launch the GUI via the command:                                                *
 *                                                                                *
 *    root -l ./TMVAGui.C                                                         *
 *                                                                                *
 **********************************************************************************/

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVAGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#endif

using namespace std;

int TMVAClassification_main(TString variable_name, TString type,int nTrees = 500, float MinNodeSize = 0.05, int maxDepth = 4)
{
	TString myMethodList = ""; 
   // The explicit loading of the shared libTMVA is done in TMVAlogon.C, defined in .rootrc
   // if you use your private .rootrc, or run from a different directory, please copy the
   // corresponding lines from .rootrc

   // methods to be processed can be given as an argument; use format:
   //
   // mylinux~> root -l TMVAClassification.C\(\"myMethod1,myMethod2,myMethod3\"\)
   //
   // if you like to use a method via the plugin mechanism, we recommend using
   //
   // mylinux~> root -l TMVAClassification.C\(\"P_myMethod\"\)
   // (an example is given for using the BDT as plugin (see below),
   // but of course the real application is when you write your own
   // method based)

   //---------------------------------------------------------------
   // This loads the library
   TMVA::Tools::Instance();

   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;

   Use["Fisher"]          = 0;
   Use["MLPBNN"]          = 0; // Recommended ANN with BFGS training method and bayesian regulator
   // --- Boosted Decision Trees
   Use["BDTB"]            = 0; // decorrelation + Adaptive Boost
   Use["BDTG"]            = 0; 
   Use["BDT"]            = 1; 
   // 
   // ---------------------------------------------------------------

   std::cout << std::endl;
   std::cout << "==> Start TMVAClassification" << std::endl;

   // Select methods (don't look at this code - not of interest)
   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

      std::vector<TString> mlist = TMVA::gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
            std::cout << std::endl;
            return 1;
         }
         Use[regMethod] = 1;
      }
   }

   // --------------------------------------------------------------------------------------------------

   // --- Here the preparation phase begins

   // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
 //  TString outfileName( "TMVA_bjet_new_powheg.root" );
 
	TString end = "axis2jet2q";
	
   TString outfileName("output/v25/TMVA_main_v25_Nm1_"+variable_name+type+end+".root" );
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
	ofstream out;
	out.open("output/v25/txt/TMVA_main_v25_Nm1_"+variable_name+type+end+".txt"); 

   // Create the factory object. Later you can choose the methods
   // whose performance you'd like to investigate. The factory is 
   // the only TMVA object you have to interact with
   //
   // The first argument is the base of the name of all the
   // weightfiles in the directory weight/
   //
   // The second argument is the output file for the training results
   // All TMVA output can be suppressed by removing the "!" (not) in
   // front of the "Silent" argument in the option string
 //  char factoryName[80];
//	sprintf(factoryName,"TMVA_13TeV_%s_%s_parameters_nTrees%i_MinNodeSize%.2f_maxDepth%i",type,end,nTrees,maxDepth);
   TMVA::Factory *factory = new TMVA::Factory( "Classification", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;P;G,D:AnalysisType=Classification" );

   // If you wish to modify default settings
   // (please check "src/Config.h" to see all available global options)
   //    (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
   
		TString weightDirName = "weights/TMVAClassification_BDTG_";
		weightDirName.Append(variable_name);
		weightDirName.Append(type);
		weightDirName.Append(end);
		weightDirName.Append("_v25");
  	   (TMVA::gConfig().GetIONames()).fWeightFileDir = weightDirName;

   // Define the input variables that shall be used for the MVA training
   // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
   // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]
  // factory->AddVariable( "myvar1 := var1+var2", 'F' );
  // factory->AddVariable( "myvar2 := var1-var2", "Expression 2", "", 'F' );



//second itereatin final order "Mqq" "DeltaPhiQQ",  "SoftN5", "HTsoft", "CSV1", "CSV2","DeltaEtaQB1", "DeltaEtaQB2",  "qgl1_VBF", "qgl2_VBF"  "Jet5_pt"

	const int max_variables_number=7;
//	TString variables_names[max_variables_number]={ "Mqq", "DeltaEtaQQ", "Jet2q_pt" ,"Jet1q_pt", "Jet1q_leadTrackPt", "Jet2q_leadTrackPt", "axis2_jet1", "axis2_jet2", "qq_pt", "Jet3_pt", "RptHard", "Zll_pt", "Zll_zstar"}; 
//	TString variables_names[max_variables_number]={ "Mqq", "DeltaEtaQQ", "Jet2q_pt" ,"Jet1q_pt", "Jet1q_leadTrackPt", "Jet2q_leadTrackPt", "axis2_jet1", "axis2_jet2", "qq_pt", "RptHard", "Zll_pt", "Zll_zstar"}; 
//	TString variables_names[max_variables_number]={ "Mqq", "DeltaEtaQQ", "Jet2q_pt",  "axis2_jet1", "axis2_jet2", "qq_pt", "RptHard", "Zll_zstar"}; 
	TString variables_names[max_variables_number]={ "Mqq", "DeltaEtaQQ","Jet2q_pt", "axis2_jet1","axis2_jet2", "qq_pt", "Zll_zstar"}; 
//	TString variables_names[max_variables_number]={ "Mqq", "DeltaEtaQQ",   "axis2_jet1" , "axis2_jet2", "Zll_zstar"}; 

 	out<<"We used variables : "<<endl;
	for (int i=0;i<max_variables_number;i++){
		if (variable_name.CompareTo("all")==0){
			factory->AddVariable(variables_names[i], "", "", 'F' );
			out<<variables_names[i]<<"   ,   ";
		}
		else {	 
			if (variables_names[i].CompareTo(variable_name)!=0) {
				factory->AddVariable(variables_names[i], "", "", 'F' );
				out<<variables_names[i]<<"   ,   ";
			}
		}
	}
	out<<endl;


   // You can add so-called "Spectator variables", which are not used in the MVA training,
   // but will appear in the final "TestTree" produced by TMVA. This TestTree will contain the
   // input variables, the response values of all trained MVAs, and the spectator variables
  // factory->AddSpectator( "spec1 := var1*2",  "Spectator 1", "units", 'F' );
  // factory->AddSpectator( "spec2 := var1*3",  "Spectator 2", "units", 'F' );

   // Read training and test data
   // (it is also possible to use ASCII format as input -> see TMVA Users Guide)
	TString fname_signal_test ="main_tmva_tree_EWK_LL_JJ_v25"+type+"_test_split_incl.root";
	TString fname_bg_test ="main_tmva_tree_DYJetstoLL_madgraph_v25"+type+"_test_split_incl.root";
	TString fname_signal_train ="main_tmva_tree_EWK_LL_JJ_v25"+type+"_train_split_incl.root";
	TString fname_bg_train ="main_tmva_tree_DYJetstoLL_madgraph_v25"+type+"_train_split_incl.root";



//	TString fname = "/afs/cern.ch/work/n/nchernya/Hbb/b_likelihood/spring15_M125_powheg_VBF_b_likelihood.root";	

   if (gSystem->AccessPathName( fname_signal_train )) { // file does not exist in local directory
		cout<<"input file "<< fname_signal_train<<" doesn't exist!"<<endl;
		//break;
	}
   if (gSystem->AccessPathName( fname_bg_train)) { // file does not exist in local directory
		cout<<"input file "<< fname_bg_train<<" doesn't exist!"<<endl;
		//break;
	}
    //  gSystem->Exec("wget http://root.cern.ch/files/tmva_class_example.root");
   
   TFile *input_signal_test = TFile::Open( fname_signal_test );
   TFile *input_bg_test = TFile::Open( fname_bg_test );
   TFile *input_signal_train = TFile::Open( fname_signal_train );
   TFile *input_bg_train = TFile::Open( fname_bg_train );
   
   std::cout << "--- TMVAClassification       : Using input signal file: " << input_signal_train->GetName() << std::endl;
   std::cout << "--- TMVAClassification       : Using input bg file: " << input_bg_train->GetName() << std::endl;
   
   // --- Register the training and test trees

   TTree *signal_train     = (TTree*)input_signal_train->Get("TMVA");
   TTree *bg_train = (TTree*)input_bg_train->Get("TMVA");
   TTree *signal_test     = (TTree*)input_signal_test->Get("TMVA");
   TTree *bg_test = (TTree*)input_bg_test->Get("TMVA");
   
   // global event weights per tree (see below for setting event-wise weights)
   Double_t signalWeight     = 1.0;
   Double_t bgWeight = 1.0;
   
   // You can add an arbitrary number of signal or background trees
   factory->AddSignalTree    ( signal_train,     signalWeight,TMVA::Types::kTraining  );
   factory->AddBackgroundTree( bg_train, bgWeight ,TMVA::Types::kTraining);
   factory->AddSignalTree    ( signal_test,     signalWeight,TMVA::Types::kTesting  );
   factory->AddBackgroundTree( bg_test, bgWeight ,TMVA::Types::kTesting);
	factory->SetSignalWeightExpression("weight");
	factory->SetBackgroundWeightExpression("weight");
   TCut mycuts = ""; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
   TCut mycutb = ""; // for example: TCut mycutb = "abs(var1)<0.5";
   


   // Tell the factory how to use the training and testing events
   //
   // If no numbers of events are given, half of the events in the tree are used 
   // for training, and the other half for testing:
   //    factory->PrepareTrainingAndTestTree( mycut, "SplitMode=random:!V" );
   // To also specify the number of testing events, use:
   //    factory->PrepareTrainingAndTestTree( mycut,
   //                                         "NSigTrain=3000:NBkgTrain=3000:NSigTest=3000:NBkgTest=3000:SplitMode=Random:!V" );
   ////////factory->PrepareTrainingAndTestTree( mycuts, mycutb,
           //                             "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );
  factory->PrepareTrainingAndTestTree( mycuts, mycutb,"!V" );

   // ---- Book MVA methods
   //
   // Please lookup the various method configuration options in the corresponding cxx files, eg:
   // src/MethoCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.html
   // it is possible to preset ranges in the option string in which the cut optimisation should be done:
   // "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable

   // Fisher discriminant (same as LD)
   if (Use["Fisher"])
      factory->BookMethod( TMVA::Types::kFisher, "Fisher", "H:!V:Fisher:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );

   if (Use["MLPBNN"])
      factory->BookMethod( TMVA::Types::kMLP, "MLPBNN", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:UseRegulator" ); // BFGS training with bayesian regulators


   if (Use["BDTD"]) // Decorrelation + Adaptive Boost
      factory->BookMethod( TMVA::Types::kBDT, "BDTD",
                           "!H:!V:NTrees=400:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:VarTransform=Decorrelate" );


	char parString[150];
	 sprintf(parString,"!H:!V:NTrees=%i:MinNodeSize=%.2f:BoostType=Grad:Shrinkage=.1:UseBaggedBoost:BaggedSampleFraction=0.6:MaxDepth=%i:nCuts=20:NegWeightTreatment=IgnoreNegWeightsInTraining",nTrees,MinNodeSize,maxDepth);
//	 sprintf(parString,"!H:!V:NTrees=%i:MinNodeSize=%.2f:MaxDepth=%i:BoostType=AdaBoost:AdaBoostBeta=0.1:SeparationType=GiniIndex:PruneMethod=NoPruning",nTrees,MinNodeSize,maxDepth);
	cout<<parString<<endl;
   if (Use["BDT"]) //
      factory->BookMethod( TMVA::Types::kBDT, "BDTG",parString );
	out<<parString<<endl;
	//	   "!H:!V:NTrees=150:MinNodeSize=0.06%:BoostType=Grad:Shrinkage=.1:nCuts=20:MaxDepth=3:NegWeightTreatment=IgnoreNegWeightsInTraining:PruneMethod=NoPruning" );
	//	   "!H:!V:NTrees=200:MinNodeSize=5%:BoostType=AdaBoost:AdaBoostBeta=0.6:nCuts=-1:MaxDepth=3:NegWeightTreatment=IgnoreNegWeightsInTraining:PruneMethod=NoPruning" );
   
	// For an example of the category classifier usage, see: TMVAClassificationCategory

   // --------------------------------------------------------------------------------------------------

   // ---- Now you can optimize the setting (configuration) of the MVAs using the set of training events

   // ---- STILL EXPERIMENTAL and only implemented for BDT's ! 
   // factory->OptimizeAllMethods("SigEffAt001","Scan");
   // factory->OptimizeAllMethods("ROCIntegral","FitGA");

   // --------------------------------------------------------------------------------------------------

   // ---- Now you can tell the factory to train, test, and evaluate the MVAs

   // Train MVAs using the set of training events
   factory->TrainAllMethods();

   // ---- Evaluate all MVAs using the set of test events
   factory->TestAllMethods();

   // ----- Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();

   // --------------------------------------------------------------

   // Save the output
   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVAClassification is done!" << std::endl;

   delete factory;

   // Launch the GUI for the root macros
   if (!gROOT->IsBatch()) TMVAGui( outfileName );

	out.close();
	return 0;

}
