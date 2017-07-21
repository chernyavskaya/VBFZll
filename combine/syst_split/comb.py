import csv
import itertools as IT

file_lists = [('results_split_systematics_el_amc.txt', 2), ('results_split_systematics_mu_amc.txt', 2), ('results_split_systematics_elmu_amc.txt',2)]

temp_data = []

for a_file in file_lists:
    file_h = open(a_file[0])
    a_list = []
    csv_reader = csv.reader(file_h, delimiter=',')
    for row in csv_reader:
       a_list.append(row[a_file[1]])           
    temp_data.append((n for n in a_list))
    file_h.close()

with open('output.data', 'w') as output_file:
    csv_writer = csv.writer(output_file, delimiter=',')
    for row in list(zip(*temp_data)):
        csv_writer.writerow(row)
