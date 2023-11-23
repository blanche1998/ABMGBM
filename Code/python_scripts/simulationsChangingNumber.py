#AUTHOR: Blanche Mongeon

#To run this project run       bash bash_script project when in the code folder

#we assume that the project was compiled 
#the script is there to change the parameters of the template file, save it as a new file and
#execute the project
#we loop through different simulations with some parameters being changed every time
#to run bash commands in python: 
#   import subprocess
#   subprocess.run(["ls"])

import subprocess
import numpy as np
#import scipy.io
import os
import sys
import shutil #this will be used to copy files by using the command shutil.copyfile(src, dst)
              #in that case we are copying the content of the file src into the file dst
import xml.etree.ElementTree as ET #to parse through and modify a xml file
from datetime import datetime
import csv
#from wetransfer import TransferApi
#from wetransferpy import WeTransfer
#from wetransfer import WeTransfer
import zipfile
import glob
import math
import operator

python_script_folder = os.getcwd()
code_folder = os.path.dirname(python_script_folder)
patients_folder = code_folder+'/patientsCSV'
errorFilePath = python_script_folder+"/simulationsChangingNumbers.txt"

#os._exit(os.EX_OK)

#for each of the patients, we simulate 6 treatment combo: No treatment, TMZ, TMZ+ICI,TMZ+OV,TMZ+ICI+OV, and ICI+OV
treatments = ["TMZ+ICI+OV"]
Nsimul = len(treatments) #number of simulations per patient
TMZsimul = ["true"]
ICIsimul = ["true"]
OVsimul = ["true"]

Nreplicate = 3 #number of replicates per simulation

#get the arguments that were passed 
arg = sys.argv
if len(arg)>1:
    executable = sys.argv[1] #this will be the project called for all simulations
else:
    #define here the names of the executables to be called for each simulation (1xN string array)
    #please consider that all projects should have been built already, and the template files will differ
    #executable = ['gbm_tmz','gbm_tmz','gbm_tmz']
    sys.exit('You must enter the name of the executable project for all simulations.')

#os.chdir(code_folder) #Go to the PhysiCell folder
xml_extension = ".xml"
template_file = "PhysiCell_settings_TMZ_ICI_OV"

os.chdir(code_folder) #make sure we're in the Code folder

patient = ["Cancer2000","Cancer3000","Cancer4000","Cancer5000","Cancer6000","Cancer7000","Cancer8000","Cancer9000","Cancer10000",
           "CTL2000","CTL3000","CTL4000","CTL5000","CTL6000","CTL7000","CTL8000","CTL9000","CTL10000",
           "TH2000","TH3000","TH4000","TH5000","TH6000","TH7000","TH8000","TH9000","TH10000",
           "Macrophages2000","Macrophages3000","Macrophages4000","Macrophages5000","Macrophages6000","Macrophages7000","Macrophages8000","Macrophages9000","Macrophages10000",
           "Stroma2000","Stroma3000","Stroma4000","Stroma5000","Stroma6000","Stroma7000","Stroma8000","Stroma9000","Stroma10000"]

for nPatient in range(len(patient)):
    px_nbr = patient[nPatient]#"topRepCell"#"CellsCSV"
    px_IC = patient[nPatient]+".csv"
    errorFile = open(errorFilePath,"a")
    now = datetime.now()
    date_string = now.strftime("%d/%m/%Y %H:%M")
    errorFile.write("\n"+date_string+": Starting simulations")
    errorFile.close()
    os.mkdir("output_"+px_nbr)
    patient_path = code_folder+"/output_"+px_nbr
    #loop through your simulations
    #i is the simulation index
    length = len(px_nbr)
    Str2 = operator.getitem(px_nbr, slice(length-4, length))
    C0 = int(Str2)

    for i in range(Nsimul):
        
        #keep track of what simulation you are starting and when
        now = datetime.now()
        date_string = now.strftime("%d/%m/%Y %H:%M")
        errorFile = open(errorFilePath,"a")
        errorFile.write("\n"+date_string+": now running simulation "+str(i+1)+" out of "+str(Nsimul)+" for "+px_nbr)
        errorFile.write('\n')
        errorFile.close()
        
        errorFile = open(errorFilePath,"a")
        #each time we start with a new type of simulation we reset the project since we will be doing many replicates
        subprocess.run(["make","clean"],stdout=errorFile)
        errorFile.write('\n')
        subprocess.run(["make","reset"],stdout=errorFile)
        errorFile.write('\n')
        subprocess.run(["make","GBM-immune-TMZ-ICI-OV"],stdout=errorFile)
        errorFile.write('\n')
        subprocess.run(["make"],stdout=errorFile)
        errorFile.write("\n Project created for patient "+str(px_nbr)+" simulation "+str(i+1))
        errorFile.close()

        shutil.copyfile('patientsCSV/'+px_IC,'cells.csv')
        shutil.copyfile('patientsCSV/'+px_IC,'config/cells.csv')

        file_name_main = treatments[i]
        os.mkdir("output_"+px_nbr+"/"+file_name_main)
        treatment_folder = code_folder+"/output_"+px_nbr+"/"+file_name_main
        #os.mkdir(treatment_folder)
        for nR in range(Nreplicate):    
            #break
            file_name = file_name_main+"_"+px_nbr+"_rep"+str(nR+1)
            output_dir = "output_"+px_nbr+"/"+file_name_main+"/"+"rep"+str(nR+1)

            #if the output directory already exists, we want to remove it
            #Be careful to make sure that you have already saved all of the output directory that
            #you want to keep
            if os.path.isdir(output_dir):
                shutil.rmtree(output_dir)

            #create the output directory
            #os.mkdir("output_"+px_nbr)
            #os.mkdir("output_"+px_nbr+"/"+file_name_main)
            os.mkdir(output_dir)

            #change the working directory
            os.chdir('config/')

            #copy the template file into the new project file
            #shutil.copyfile('~/PhysiCell/sample_projects/'+proj+'/config'+template_file+xml_extension,file_name+xml_extension)
            shutil.copyfile(template_file+xml_extension,file_name+xml_extension)

            #-------------this is where you parse through and modify your xml file as needed-------------
            tree = ET.parse(file_name+xml_extension)
            root = tree.getroot()
            rootTag = root.tag

            #set the maximal simulation time 

            userparam_element = root.find('user_parameters')
            TMZ_tx = userparam_element.find("TMZ_tx")
            TMZ_tx.text = TMZsimul[i]
            ICI_tx = userparam_element.find("ICI_tx")
            ICI_tx.text = ICIsimul[i]
            OV_tx = userparam_element.find("OV_tx")
            OV_tx.text = OVsimul[i]

            V0 = userparam_element.find("V0")

            #to calculate the new V0, we need to know how many cells this patient has in its sample
            #this loop is run for 
            newV0 = ((5/9740)*C0*C0)/(20*20*math.pi)

            V0.text = str(newV0)

            #change the output file
            save_element = root.find('save')
            save_element.find('folder').text = output_dir

            #save the resulting xml file
            tree.write(file_name+xml_extension)

            #---------------------------------------Run the project--------------------------------------

            #go back to the PhysiCell directory
            #the console output of this specific simulation will be saved to a specific file
            os.chdir(code_folder)
            outfile_name = 'outfile_'+file_name
            with open(outfile_name, "w") as outfile:
                subprocess.run(["./"+executable,"config/"+file_name+xml_extension],stdout=outfile)

    errorFile = open(errorFilePath,"a")
    errorFile.write("\n "+px_nbr+" completed!")
    errorFile.close()
