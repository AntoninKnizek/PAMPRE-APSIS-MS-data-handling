import numpy as np
import os
import database_functions as df
import shutil

# range of m/z in database 
mzmin = 1 #minimum m/z in the data
mzmax = 200 #maximum m/z in the data
dm = mzmax-mzmin+1 # number of points in the data
mz = np.linspace(mzmin,mzmax,dm,dtype=int) # array of m/z for the data
NIST_files = './NIST'

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

db = df.add_compound(db,helium,He,0.326,4,mzmin,mzmax)

db = df.add_fragment(db,helium,4,100)

db = df.add_compound_from_NIST(db,'./NIST/74-82-8-Mass.jdx')
'''The file used in the above function can be downloaded from NIST chemistry 
webbook - methane.'''


db = df.normalize_intensity(db) #normalizes intensities in the database

print('This is the current database.')
print(db)

############################################################
shutil.copy('central_database.txt', 'central_database_2.txt')
df.add_compound_from_NIST_to_file('central_database_2.txt',
                                        './NIST/74-82-8-Mass.jdx')
                                        
print('A new compound was added to the central database.')

###########################################################
compounds = df.identify_by_mz(db, 25)
print(*compounds)
