'''This code is written to prepare database for MS analysis.'''

import numpy as np
import os
import database_functions as df

##### BEGINNING OF THE PART TO BE MODIFIED  

# range of m/z in database 
mzmin = 1 #minimum m/z in the data
mzmax = 200 #maximum m/z in the data
dm = mzmax-mzmin+1 # number of points in the data
mz = np.linspace(mzmin,mzmax,dm,dtype=int) # array of m/z for the data

output_path = '.' # path directory where the database text file will be created 

##### END OF THE PART TO BE MODIFIED 

################## open file database.txt
db = df.read_central_database('central_database.txt',mzmin,mzmax)
db = df.choose_from_central_database(db,*('hydrogen',  
                                          'water', 
                                          'argon', 
                                          'methane', 
                                          'acetylene', 
                                          'ethylene', 
                                          'ethane', 
                                          'allene', 
                                          'benzene', 
                                          'propane', 
                                          'diacetylene', 
                                          'propyne', 
                                          '1,3-butadiene', 
                                          '2-methylpropane', 
                                          '2-methylbutane', 
                                          'propene', 
                                          'toluene', 
                                          'isohexane', 
                                          'n-butane', 
                                          'cis-2-butene', 
                                          '1-butene', 
                                          'trans-2-butene', 
                                          'o-xylene', 
                                          '3-methyl-1-butene', 
                                          'hex-1-ene', 
                                          '2.2-dimethylpropane', 
                                          'trans-1,3-pentadiene', 
                                          '1,3,5-heptatriene', 
                                          '2,3-dimethylpent-2-ene'
                                          )
                                      )

db = df.normalize_intensity(db) #normalizes intensities in the database

# remove the file if it already exists
if os.path.isfile(f'{output_path}/database.txt') == True:
    os.remove(f'{output_path}/database.txt')
 
database = open(f'{output_path}/database.txt', 'a')

# first line of the file 
database.write('Compound_index'+'\t')
database.write('Compound_name'+'\t')
database.write('Compound_Formula'+'\t')
database.write('ionization_cross_section'+'\t')
database.write('mz_simple_ioni'+'\t')
for ti in range(dm):
    if ti<dm-1:
        database.write('mz'+str(mz[ti])+'\t')
    else:
        database.write('mz'+str(mz[ti]))

# other lines with parameters and database mass spectrum of each gas 
for i,key in enumerate(db.keys()):
    database.write('\n')
    database.write(str(i+1)+'\t')
    database.write(str(db[key]['name'])+'\t')
    database.write(str(db[key]['formula'])+'\t')
    database.write(str(db[key]['ioniz_cross_sec'])+'\t')
    database.write(str(int(db[key]['M']))+'\t')
    for a in range(dm):
        if a<dm-1:
            database.write(str(db[key]['intensity'][a])+'\t')
        else:
            database.write(str(db[key]['intensity'][a]))

###########################################
# end of the file 
database.close() # close file 
