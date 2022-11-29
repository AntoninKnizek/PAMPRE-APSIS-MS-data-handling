import pylab as plt
import data_handling as dh

filename = '2022_10_26_bargraph_gases_after_experiment.csv'

data_csv = dh.load_bargraph_csv(f'./data/{filename}')

data_dict = dh.get_bargraph_dict_from_array(data_csv)

data_avg = dh.get_bargraph_average(data_dict)

plt.ion()
plt.title(f'{filename}')
plt.xlabel('m/z')
plt.ylabel('Intensity (a.u.)')
plt.yscale('log')

plt.bar(data_avg[:,0],data_avg[:,1])
plt.show()
