from itertools import izip # like zip but gives us an iterator
import sys

dy_opt = sys.argv[1]

with open('results_split_systematics_el_%s.txt'%dy_opt) as f1, open('results_split_systematics_mu_%s.txt'%dy_opt) as f2,open('results_split_systematics_elmu_%s.txt'%dy_opt) as f3, open('unc_splitting_%s_exel.txt'%dy_opt, 'w') as out:
    line_counter=0
    for f1line, f2line, f3line in izip(f1, f2, f3):
        if line_counter==0 : 
          ch = ','
        else :
          ch = '+'
        f1line_start=f1line[:f1line.find(',')]
        f1line_end=f1line[f1line.find(ch):]
        f2line=f2line[f2line.find(ch):]
        f3line=f3line[f3line.find(ch):]
         # out.write('{}\t\t\t{}\t\t\t{}'.format(f1line.strip(), f2line.strip(),f3line))
        out.write('%s\t\t\t\t %s\t\t\t\t %s\t\t\t\t %s\n'%(f1line_start,f1line_end.strip(), f2line.strip(), f3line.strip()))
        line_counter+=1
