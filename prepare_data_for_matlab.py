import pylab as plt
import data_handling as dh

filename = '2022_10_26_bargraph_gases_after_experiment.csv'

data_csv = dh.load_bargraph_csv(f'./data/{filename}')

data_dict = dh.get_bargraph_dict_from_array(data_csv)

data_avg = dh.get_bargraph_average(data_dict)

data_avg = dh.normalize_data(data=data_avg,mz=26)

with open('data_bar.txt', 'w') as f:
    for i in range(len(data_avg[:,1])):
        f.write(str(data_avg[i,1]))
        f.write('\n')
f.close()
