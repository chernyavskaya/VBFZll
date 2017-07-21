import os
import sys
type = sys.argv[1]
dy_option = sys.argv[2]
intOption="_interference"
ntot=0
types=['mu','el','elmu']

output_txt="results_split_systematics_%s_%s.txt"%(type,dy_option)
output = open(output_txt,"wt")
#output.write('\t\t\t\t\t\t\t\t\t\t\tmu\t\t\t\t\t\t\t\t\t\t\t\tel\t\t\t\t\t\t\t\t\t\t\tel and mu')
output.write('Unc , %s\n'%type)
if True:
#for type in types : 
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
# for     	
   for current_syst in range(-2,ntot):   #loop over systematics
 #  for current_syst in range(-2,3):   #loop over systematics
          card_name="../cards/ewkZjj_13TeV_datacard"+ intOption + "_" + type + "_DY"+dy_option+"_atanh_v25alldata.txt"
          card = open(card_name, "rt")
          lines_card=card.readlines()
          input_txt = 'exp_fit_'+type+"_DY"+dy_option+"_v25alldata_%i.txt"%current_syst
          result = open(input_txt, "rt")
          print 'current_syst',current_syst
          lines_result=result.readlines()
          line_counter=0
          if current_syst==-2 :
            text = 'syst+stat'
            text2 =lines_result[0]
            text2 =text2[text2.find('-'):text2.find('(')]
           # output.write('%-35s %-35s\n'%(text,text2))
            output.write('%s,%s\n'%(text,text2))
            continue
          if current_syst==-1 :
            text = 'stat'
            text2 =lines_result[0]
            text2 =text2[text2.find('-'):text2.find('(')]
          #  output.write('%-35s %-35s\n'%(text,text2))
            output.write('%s,%s\n'%(text,text2))
            continue
          else:
            text = lines_card[s_start+current_syst]
            if type!='elmu' : text = text[:text.find('\t')]
            else : text = text[:text.find(' ')]
            text2 =lines_result[0]
            text2 =text2[text2.find('-'):text2.find('(')]
         #   output.write('%-35s %-35s\n'%(text,text2))
            output.write('%s,%s\n'%(text,text2))



