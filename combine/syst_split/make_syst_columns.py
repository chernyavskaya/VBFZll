from itertools import izip # like zip but gives us an iterator

with open('results_split_systematics_el_amc.txt') as f1, open('results_split_systematics_mu_amc.txt') as f2,open('results_split_systematics_elmu_amc.txt') as f3, open('output.txt', 'w') as out:
    for f1line, f2line, f3line in izip(f1, f2, f3):
        f1line_start=f1line[:f1line.find(',')]
        f1line_end=f1line[f1line.find('-'):]
        f2line=f2line[f2line.find('-'):]
        f3line=f3line[f3line.find('-'):]
       # out.write('{}\t\t\t{}\t\t\t{}'.format(f1line.strip(), f2line.strip(),f3line))
        out.write('%-30s %-30s %-30s %-30s\n'%(f1line_start,f1line_end.strip(), f2line.strip(), f3line.strip()))
