import os
import sys

combine_line='combine -M MaxLikelihoodFit -t -1 --expectSignal=1 --robustFit=1 --stepSize=0.05 --preFitValue=1.0 --rMin=-5 --rMax=5'
combine_line_nosyst='combine -M MaxLikelihoodFit  --X-rtd FITTER_NEW_CROSSING_ALGO  -t -1 --expectSignal=1 --robustFit=1 --preFitValue=1.0  --stepSize=0.05 --rMin=-5 --rMax=5 -S 0'

type = sys.argv[1]
dy_option = sys.argv[2]
intOption="_interference"
ntot=0
if type!='elmu' : 
	s_start=14
	if dy_option=='amc' :
		ntot = 15 
	if dy_option=='mdg' :
		ntot = 16
else :
	s_start=16
	if dy_option=='amc' :
		ntot = 16 
	if dy_option=='mdg' :
		ntot = 17
	
if (1>0) : 
   for current_syst in range(-2,ntot):   #loop over systematics
      card_name="ewkZjj_13TeV_datacard"+ intOption + "_" + type + "_DY"+dy_option+"_atanh_v25alldata_%i.txt"%current_syst
      with open(card_name, "wt") as card:
        output_txt = 'exp_fit_'+type+"_DY"+dy_option+"_v25alldata_%i.txt"%current_syst
        template_name="../cards/ewkZjj_13TeV_datacard"+ intOption + "_" + type + "_DY"+dy_option+"_atanh_v25alldata.txt"
        with open(template_name, "rt") as template:
          line_counter=0
          if current_syst==-2 :
               print 'current systematic', current_syst
               print combine_line,"  ", template_name   
               os.system("%s %s  | grep 'Best fit r: 1' > %s  "%(combine_line,template_name,output_txt))   #run combine and save output in some txt file
               continue
          if current_syst==-1 :
               print 'current systematic', current_syst
               print combine_line_nosyst,"  ", template_name   
               os.system("%s %s  | grep 'Best fit r: 1' > %s  "%(combine_line_nosyst,template_name,output_txt))   #run combine and save output in some txt file
               continue
          for line in template:
            if line_counter < s_start :
              card.write(line)
              line_counter+=1
              continue
            else:
              line_counter+=1
              if ('_stats_DY' in line ) and (current_syst==ntot-1) :
                line = '#' + line
              elif (line_counter-current_syst-1 == s_start) :
                line = '#' + line
              card.write(line)
          card.close()
          print 'current systematic', current_syst
          os.system("%s %s  | grep 'Best fit r: 1' > %s  "%(combine_line,card_name,output_txt))   #run combine and save output in some txt file



