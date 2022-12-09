Python MS data analysis library for PAMPRE/APSIS
Author: Antonín Knížek (Antonin Knizek), email: knizeka@gmail.com
Created 29 Nov 2022.

Introduction:
This library was craeted to facilitate mass spectrometry data handling and analaysis for data generated on PAMPRE and APSIS experiments at LATMOS, UVSQ, Paris, France.
The library can in principle work for any data and for data generated on Hiden Analytical spectrometers especially.
This intended use of this library is preparation of data inputs for the matlab code for MS data analysis developped by Thomas Gautier (DOI: 10.1002/rcm.8684).

#######################################################
This library now contains the following files:
- data_handling.py
- database_functions.py
- central_database.py
- prepare_data_for_matlab.py
- plot_bargraph.py
- create_database_for_matlab.py
- README.txt
- LICENSE

#######################################################
Prerequisites for use:
    This library was written on Ubuntu 22.04 LTS with python 3.10.6 and was not tested on other versions or operating systems
    Required packages: numpy, json, random, matplotlib

#######################################################
Description of the files:
data_handling.py:
    This module is intended for loading and manipulating the experimental data.

    Functions:
        load_bargraph_csv(path_to_file): 
            This function loads an experimental csv file and returns a data array.
        
        get_bargraph_dict_from_array(data_arr,start_cycle,end_cycle):
            This function takes in a data array generated from load_bargraph_csv or one in identical format and products a dictionary with meaasurement cycles as the dictionary keys
            start_cycle and end_cycle are optional parameters for specyfing only some measurement cycles to use.
        
        get_bargraph_average(data_dict):
            This function averages the experimental data across the cycles in the dictionary from get_bargraph_dict_from_array
        
        normalize_data(data,mz):
            This function normalizes the experimental data in the form of a 2D array (from the get_bargraph_average) to a specified mz value.
            
            
prepare_data_for_matlab.py
    A short example script that uses functions from the data_handling.py module and produces data_bar.txt file, which can be directly imported to the matlab code.
    
plot_bargraph.py
    A short example script ths uses functions from the data_handling.py module and creates a spectrum plot in matplotlib.

central_database.txt:
    This is a central database with molecules, their metadata and fragmentation patterns.
    This database is in json-compatible format.

database_functions.py:
    This module is intended for loading the central database/creating a database for use for the matlab code.
    This database contains functions to read a central database stored as central_database.txt 

    Functions:
        add_compound(db,name,form,ics,M,mz_min,mz_max):
            This functions adds a molecule with all its metadata provided as arguments to a dictionary database (argument db).
            The arguments are:
                db - name of the database
                name - compound name
                formula - compound chemical formula
                ics - ionization cross section in square Ångström
                M - molar mass in g mol-1
                int_min - minimum mz
                int_max - maximum mz
                
        add_fragment(db,compound,mz,i):
            Function that adds a fragment intensity to the database in live code.
            db - name of the database
            compound - compound name
            mz - m/z of the fragment, format = integer
            i - fragment intensity in counts
                
        normalize_intensity(db):
            Function to normalize the intensities in the database.
            db - database name
            
        read_central_database(database_name,mz_min,mz_max):
            This function loads the central database and adapts the intensity field
            to correct length.
            database name = filename of the database
            mz_min - min m/z for the experiment
            mz_max - max m/z for the experiment
            
        choose_from_central_database(c_db,*compounds):
            This function extracts a subset of molecules from the central database.
            c_db - central database
            db - the new database
            compounds - tuple with compound names to be extracted.
            
        identify_by_mz(c_db, mz):
            This function lists molecules from the central database
             with fragments at the required mass mz.
            c_db - central database
            mz - required mz
            
        add_compound_from_NIST(db,path_to_file):
            This function will load a NIST jcampdx mass spectrum and add it to the 
            live database.
            db - database
            path_to_file : path to NIST jdx file
            
        add_compound_from_NIST_to_file(path_to_c_db,path_to_compound):
            This function will read a NIST jdx file and write the compound
            data to the central database file.
            path_to_c_db - path to the central database file
            path_to_comopund - path to the jdx file
            WARNING: This function overwrites the central database file. Keep a backup
            copy, save lives.
            WARNING: This function does not check for duplicates in the database.
        round_half_up(n,decimals):
            a funciton for rounding up numbers
            
create_database_for_matlab.py:
    A short example code which loads the central database, normalizes the fragmentation pattern intesities and creates a database.txt file, which can be directly improted by the matlab code.
    Also creates a constraint file with data from an IR analysis to which this code can be coupled. For this purpose, modify the relative concentrations and add/remove compounds.
    
    
example_database.py:
    A short example code which shows how all functions from the database_functions.py module are used.
    
#######################################################           
Credits:
This library was created by Antonín Knížek.
The matlab code for which the database is intended was written by Thomas Gautier.
The idea of the database is built on the original codes of Thomas Drant.

#######################################################
License:
See the license file.
