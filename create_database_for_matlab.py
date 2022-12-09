'''This code is written to prepare database for MS analysis.'''

import numpy as np
import os
import database_functions as df

##### BEGINNING OF THE PART TO BE MODIFIED  

# range of m/z in database 
mzmin = 1 #minimum m/z in the data
mzmax = 100 #maximum m/z in the data
dm = mzmax-mzmin+1 # number of points in the data
mz = np.linspace(mzmin,mzmax,dm,dtype=int) # array of m/z for the data

output_path = '.' # path directory where the database text file will be created 

##### END OF THE PART TO BE MODIFIED 

################## open file database.txt
db = df.read_central_database('central_database.txt',mzmin,mzmax)
db = df.choose_from_central_database(db,*('water',
                                          'carbon dioxide',
                                          'argon',
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
                                          'cis-2-butene', 
                                          '1-butene', 
                                          'trans-2-butene', 
                                          '3-methyl-1-butene', 
                                          'hex-1-ene', 
                                          'trans-1,3-pentadiene', 
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
    
    r = int(df.round_half_up(float(db[key]['M'])))
    database.write(str(r)+'\t')
    for a in range(dm):
        if a<dm-1:
            database.write(str(db[key]['intensity'][a])+'\t')
        else:
            database.write(str(db[key]['intensity'][a]))

###########################################
# end of the file 
database.close() # close file 


###########################################
#OPTIONAL: print out relative flow factors in temrinal:
flows = []

for key in db.keys():
    db[key]['flow'] = np.sqrt(db['water']['M']/db[key]['M'])
    flows.append(db[key]['flow'])
print(f'flow = {flows}')

###########################################
#OPTIONAL: print a gain vector (unknown gains set to 1):
gains = []
for key in db.keys():
    if db[key]['gain'] == '':
        db[key]['gain'] = 1
    gains.append(db[key]['gain'])
print(f'gain = {gains}')

############################################
#OPTIONAL: create file with the constraints for the lsqlin solver:

#data from IR:
#c(acetylene) = 1
#c(acetylene)/c(ethane) = 136.894
#c(acetylene)/c(ethylene) = 152.312
#c(acetylene)/c(CO2) = 31.112
#c(acetylene)/c(water) = 2.098

###Modify the database:
db['acetylene']['intensity'] = (
     db['acetylene']['intensity']
     
    + (db['acetylene']['intensity']
      *db['acetylene']['gain']
      *db['acetylene']['flow']
      *db['acetylene']['ioniz_cross_sec']
      /db['ethane']['gain']
      /db['ethane']['flow']
      /db['ethane']['ioniz_cross_sec']
      /136.894
      )
      
    + (db['acetylene']['intensity']
      *db['acetylene']['gain']
      *db['acetylene']['flow']
      *db['acetylene']['ioniz_cross_sec']
      /db['ethylene']['gain']
      /db['ethylene']['flow']
      /db['ethylene']['ioniz_cross_sec']
      /152.312
      )
      
    + (db['acetylene']['intensity']
      *db['acetylene']['gain']
      *db['acetylene']['flow']
      *db['acetylene']['ioniz_cross_sec']
      /db['carbon dioxide']['gain']
      /db['carbon dioxide']['flow']
      /db['carbon dioxide']['ioniz_cross_sec']
      /31.112
      )
      
    + (db['acetylene']['intensity']
      *db['acetylene']['gain']
      *db['acetylene']['flow']
      *db['acetylene']['ioniz_cross_sec']
      /db['water']['gain']
      /db['water']['flow']
      /db['water']['ioniz_cross_sec']
      /2.098
      )
      
  )


db['ethane']['intensity'] = np.zeros(dm)
db['ethylene']['intensity'] = np.zeros(dm)
db['carbon dioxide']['intensity'] = np.zeros(dm)
db['water']['intensity'] = np.zeros(dm)

####Write the constraints matrix file:
# remove the file if it already exists
if os.path.isfile(f'{output_path}/constraint.txt') == True:
    os.remove(f'{output_path}/constraint.txt')
 
constraint = open(f'{output_path}/constraint.txt', 'a')

# first line of the file 
constraint.write('Compound_index'+'\t')
constraint.write('Compound_name'+'\t')
constraint.write('Compound_Formula'+'\t')
constraint.write('ionization_cross_section'+'\t')
constraint.write('mz_simple_ioni'+'\t')
for ti in range(dm):
    if ti<dm-1:
        constraint.write('mz'+str(mz[ti])+'\t')
    else:
        constraint.write('mz'+str(mz[ti]))

# other lines with parameters and database mass spectrum of each gas 
for i,key in enumerate(db.keys()):
    constraint.write('\n')
    constraint.write(str(i+1)+'\t')
    constraint.write(str(db[key]['name'])+'\t')
    constraint.write(str(db[key]['formula'])+'\t')
    constraint.write(str(db[key]['ioniz_cross_sec'])+'\t')
    
    r = int(df.round_half_up(float(db[key]['M'])))
    constraint.write(str(r)+'\t')
    for a in range(dm):
        if a<dm-1:
            constraint.write(str(db[key]['intensity'][a])+'\t')
        else:
            constraint.write(str(db[key]['intensity'][a]))

###########################################
# end of the file 
constraint.close() # close file 


