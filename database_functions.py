import numpy as np
import json
import random
import math

def add_compound(db,name,form,ics,M,mz_min,mz_max,gain=None,flow=None,
                 source=None, doi=None,sensitivity=None,comments=None):
    '''Function to add new species to a database in live code.
    db - name of the database
    name - compound name
    formula - compound chemical formula
    ics - ionization cross section in square Ångström
    M - molar mass in g mol-1
    int_min - minimum mz
    int_max - maximum mz
    '''
    if name not in db.keys():
        db[name] = {}
        db[name]['name'] = name
        db[name]['formula'] = form
        db[name]['ioniz_cross_sec'] = ics
        db[name]['M'] = M
        db[name]['m/z'] = [*range(mz_min,mz_max+1)]
        db[name]['intensity'] = np.zeros(len(db[name]['m/z']), dtype=int)
        db[name]['gain'] = gain
        db[name]['flow'] = flow
        db[name]['source'] = source
        db[name]['doi'] = doi
        db[name]['sensitivity'] = sensitivity
        db[name]['comments'] = comments
    else:
        print('This compound already exists in the database.')
    
    return db
    
def add_fragment(db,compound,mz,i):
    '''Function that adds a fragment intensity to the database in live code.
    db - name of the database
    compound - compound name
    mz - m/z of the fragment, format = integer
    i - fragment intensity in counts
    '''
    if compound in db.keys():
        pos_index = db[compound]['m/z'].index(mz)
        db[compound]['intensity'][pos_index] = i
    
    else:
        'This compound does not exist in the database.'
    
    return db

def normalize_intensity(db):
    '''Function to normalize the intensities in the database.
    db - database name
    '''
    
    for i in db.keys():
        db[i]['intensity'] = ((db[i]['intensity'])
                              /np.max(db[i]['intensity'])
                              *9999
                             )
    return db
    
def read_central_database(database_name,mz_min,mz_max):
    '''This function loads the central database and adapts the intensity field
    to correct length.
    database name = filename of the database
    mz_min - min m/z for the experiment
    mz_max - max m/z for the experiment
    '''
    
    f = open(f'./{database_name}','r')
    c_db = json.load(f)
    f.close()
    
    for key in c_db.keys(): #loops through compounds
        c_db[key]['m/z'] = [*range(mz_min,mz_max+1)]
        a = np.zeros(len(c_db[key]['m/z']))
        for val in c_db[key]['intensity'].keys(): #loops through given fragments
            if int(val) in c_db[key]['m/z']:
                a[c_db[key]['m/z'].index(int(val))] = (
                                                    c_db[key]['intensity'][val])
        c_db[key]['intensity'] = list(a)
    
    return c_db
    
def choose_from_central_database(c_db,*compounds):
    '''This function extracts a subset of molecules from the central database.
    c_db - central database
    db - the new database
    compounds - tuple with compound names to be extracted.'''
    
    db = {}
    for compound in compounds:
        if compound not in c_db.keys():
            print(f'{compound} is not in the database.')
        else:
            db[compound] = c_db[compound]
            
    return db
    
def identify_by_mz(c_db, mz):
    '''This function lists molecules from the central database
     with fragments at the required mass mz.
    c_db - central database
    mz - required mz
    '''
    
    compounds = []
    for compound in c_db.keys():
        i = c_db[compound]['m/z'].index(mz)
        
        if c_db[compound]['intensity'][i] != 0:
            compounds.append(c_db[compound]['name'])
            
    return tuple(compounds)
    
def add_compound_from_NIST(db,path_to_file):
    '''This function will load a NIST jcampdx mass spectrum and add it to the 
    live database.
    db - database
    path_to_file : path to NIST jdx file
    '''
    
    f = open(path_to_file,'r')
    lines = f.readlines()

    for ind,line in enumerate(lines):
        if 'TITLE' in line:
            title_index = ind
        if 'MOLFORM' in line:
            molform_index = ind
        if 'MW' in line:
            mw_index = ind
        if 'PEAK TABLE' in line:
            pks_start = ind+1
        if 'END=' in line:
            pks_end = ind-1
            
    pks_range = range(pks_start,pks_end)
    
    '''This block reads desired parameters except the spectrum.'''
    name = lines[title_index][lines[title_index].rfind('=')+1:-1]
    M = lines[mw_index][lines[mw_index].rfind('=')+1:-1]
    molform = lines[molform_index][lines[molform_index].rfind('=')+1:-1]
    molform = molform.replace(' ','')
    
    '''This block, however clumsily, reads the spectrum itself.'''    
    fragments = ''
    for num in pks_range:
        fragments = fragments+str(lines[num])
    
    fragments.replace('\n',' ')

    pairs = fragments.split(' ')
    
    mz = []
    intens = []
    
    for pair in pairs:
        mz.append(pair[:pair.index(',')])
        intens.append(pair[pair.index(',')+1:])
    
    for pos,val in enumerate(intens):
        intens[pos] = intens[pos].strip('\n')
        
    '''This block adds the compound to the live database.'''
    
    rand_key = random.choice(list(db))
    
    
    if name not in db.keys():
        db[name] = {}
        db[name]['name'] = name
        db[name]['formula'] = molform
        db[name]['ioniz_cross_sec'] = ''
        db[name]['M'] = M
        db[name]['gain'] = ''
        db[name]['flow'] = ''
        db[name]['source'] = ''
        db[name]['doi'] = doi
        db[name]['sensitivity'] = ''
        db[name]['comments'] = ''
        
        db[name]['m/z'] = db[rand_key]['m/z']
        db[name]['intensity'] = np.zeros(len(db[name]['m/z']), dtype=int)
        
        for i,ii in zip(mz,intens):
            pos_index = db[name]['m/z'].index(int(i))
            db[name]['intensity'][pos_index] = ii
    else:
        print(f'The compound {name} already exists in the database.')
        
    f.close()
    
    return db
    
    
def add_compound_from_NIST_to_file(path_to_c_db,path_to_compound):
    '''This function will read a NIST jdx file and write the compound
    data to the central database file.
    path_to_c_db - path to the central database file
    path_to_comopund - path to the jdx file
    WARNING: This function overwrites the central database file. Keep a backup
    copy, save lives.
    WARNING: This function does not check for duplicates in the database.
    '''

    f = open(path_to_compound,'r')
    lines = f.readlines()

    for ind,line in enumerate(lines):
        if 'TITLE' in line:
            title_index = ind
        if 'MOLFORM' in line:
            molform_index = ind
        if 'MW' in line:
            mw_index = ind
        if 'PEAK TABLE' in line:
            pks_start = ind+1
        if 'END=' in line:
            pks_end = ind-1
            
    pks_range = range(pks_start,pks_end)
    
    '''This block reads desired parameters except the spectrum.'''
    name = lines[title_index][lines[title_index].rfind('=')+1:-1]
    M = lines[mw_index][lines[mw_index].rfind('=')+1:-1]
    molform = lines[molform_index][lines[molform_index].rfind('=')+1:-1]
    molform = molform.replace(' ','')
    
    '''This block, however clumsily, reads the spectrum itself.'''    
    fragments = ''
    for num in pks_range:
        fragments = fragments+str(lines[num])
    
    fragments.replace('\n',' ')

    pairs = fragments.split(' ')
    
    mz = []
    intens = []
    
    for pair in pairs:
        mz.append(pair[:pair.index(',')])
        intens.append(pair[pair.index(',')+1:])
    
    for pos,val in enumerate(intens):
        intens[pos] = intens[pos].strip('\n')
    
    f.close()

    '''This block writes the spectrum to file.'''
    len_mz = len(mz)
    cb = '{'
    cb2 = '}'
    
    fin = open(path_to_c_db,'r')
    fin_lines = fin.readlines()
    fin.close()
    
    
    if '}' not in fin_lines[-1]:
        print('''Central database probably contains trailing whitespace
        at the end. Make sure the central database file end with a line 
        containing the '}' character and no more blank lines are present.''')
        
    else:
        fout = open(path_to_c_db,'w')
        for l in fin_lines[:-2]:
            fout.write(l)
        fout.write(f'   {cb2},\n')
        fout.write('\n')
        

        
        fout.write(f'  \"{name}\" : {cb}\n')
        fout.write(f'    \"name\" : \"{name}\",\n')
        fout.write(f'    \"formula\" : \"{molform}\",\n')
        fout.write(f'    \"ioniz_cross_sec\" : "",\n')
        fout.write(f'    \"M\" : {M},\n')
        fout.write(f'    \"gain\" : "",\n')
        fout.write(f'    \"flow\" : "",\n')
        fout.write(f'    \"source\" : "",\n')
        fout.write(f'    \"doi\" : "",\n')
        fout.write(f'    \"sensitivity\" : "",\n')
        fout.write(f'    \"comments\" : "",\n')
        fout.write(f'    \"intensity\" : {cb}\n')
        for a in range(len_mz-1):
            fout.write(f'        \"{mz[a]}\" : {intens[a]},\n')
        fout.write(f'        \"{mz[len_mz-1]}\" : {intens[len_mz-1]}\n')
        fout.write(f'        {cb2}\n')
        fout.write(f'    {cb2}\n')
        fout.write(f'{cb2}')
        
        fout.close()
    return print('Compound successfully added to the central database.')
    
def round_half_up(n,decimals=0):
    multiplier = 10 ** decimals
    return math.floor(n*multiplier + 0.5) / multiplier
    
