import numpy as np
import random

def load_bargraph_csv(path_to_file):
    '''Loads bargraph csv data from MS at PAMPRE into an array.'''
    data = {}
    data_list = []
    
    file = open(path_to_file, 'r')
    lines = file.readlines()
    no_lines = len(lines)
    
    for l,line in enumerate(lines):
        if l >= 42:
            line_content = line.split(';')
            for i,ele in enumerate(line_content):
                line_content[i] = ele.replace(',', '.')
                line_content[i] = ele.strip()

            data_list.append(line_content)
    no_elems = len(data_list)
    data_arr = np.array(data_list)
    file.close()
    
    return data_arr
    
def get_bargraph_dict_from_array(data_arr,start_cycle=None,end_cycle=None):
    '''This function creates a data dictionary to from the loaded csv from 
    load_bargraph_csv function.
    cycles - amount of cycles required
    '''
    
    cycles = set(data_arr[:,0])
    
    if start_cycle == None and end_cycle == None:
        cycles = cycles
    elif start_cycle == None and end_cycle != None:
        for e in cycles.copy():
            if int(e) > end_cycle:
                cycles.remove(e)
    elif start_cycle != None and end_cycle == None:
        for e in cycles.copy():
            if int(e) < start_cycle:
                cycles.remove(e)
    elif start_cycle != None and end_cycle != None:
        for e in cycles.copy():
            if int(e) > end_cycle or int(e) < start_cycle:
                cycles.remove(e)
    else:
        print('''Something went wrong with 
        the get_bargraph_dict_from_array function.''')
        
        
    data_dict = {}
    
    for c in cycles:
        data_dict[int(c)] = []

    for i in range(len(data_arr)):
        if int(data_arr[i,0]) in data_dict.keys():
            data_dict[int(data_arr[i,0])].append(data_arr[i,3:])
        
    for key in data_dict.keys():
        data_dict[key] = np.array(data_dict[key])
        
    return data_dict
    
def get_bargraph_average(data_dict):
    '''Calculates an average of the bragraph measurements 
    to give one final spectrum.'''
    avg = []
    
    ii = random.choice(list(data_dict.keys()))
    
    for mz in range(len(data_dict[ii][:,0])):
        tot = 0
        for key in data_dict.keys():
            tot += float(data_dict[key][mz,1])
        one_avg = tot/len(data_dict.keys())
        avg.append(one_avg)

    ms = []
    for m in data_dict[ii][:,0]:
        ms.append(float(m.replace(',','.')))
    
    res = np.column_stack((ms,avg))
    
    return res
    
def normalize_data(data,mz):
    '''This function normalizes data to given value.
    data - 2D data array with mz and intensities from get_bragraph_average.
    mz - value to be normalized to.
    '''
    
    data[:,1] = data[:,1]/data[list(data[:,0]).index(mz),1]
    
    return data    
