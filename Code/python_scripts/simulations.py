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

python_script_folder = os.getcwd()
code_folder = os.path.dirname(python_script_folder)
patients_folder = code_folder+'/patientsCSV'
errorFilePath = python_script_folder+"/simulations.txt"


#load the file containing the name of the patients - the file is in the python script folder 
#which we are in right now
patients = []
with open('patients_simul.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
        l = len(row)
        if line_count==0:
            for nl in range(l):
                patients.append(row[nl])
        line_count += 1

Npatient = len(patients)

errorFile = open(errorFilePath,"a")
errorFile.write("\n Patients we will simulate are:")
for npx in range(Npatient):
    errorFile.write(patients[npx]+", ")
errorFile.write('\n')
errorFile.close()

#os._exit(os.EX_OK)

#for each of the patients, we simulate 6 treatment combo: No treatment, TMZ, TMZ+ICI,TMZ+OV,TMZ+ICI+OV, and ICI+OV
treatments = ["none"]#"None","TMZ","ICI","OV","TMZ+ICI+OV","ICI+OV"
Nsimul = len(treatments) #number of simulations per patient
TMZsimul = ["false"] #"false","true","false","false","true","false"
ICIsimul = ["false"] #"false","false","true","false","true","true"
OVsimul = ["false"] #"false","false","false","true","true","true"

Nreplicate = 1 #number of replicates per simulation

#Nsimul = 3
#treatments = ["TMZ+OV","TMZ+ICI+OV","ICI+OV"]
#TMZsimul = ["true","true","false"]
#ICIsimul = ["false","true","true"]
#OVsimul = ["true","true","true"]

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

#j ios the patient's index
for j in range(0,6):
    os.chdir(code_folder) #make sure we're in the Code folder
    px_IC = patients[j] #get the patient's filename for initial conditions
    #get the patient's number
    #px_nbr = px_IC.replace('patient','')
    px_nbr = px_IC.replace('.csv','')
    
    errorFile = open(errorFilePath,"a")
    now = datetime.now()
    date_string = now.strftime("%d/%m/%Y %H:%M")
    errorFile.write("\n"+date_string+": Starting "+px_nbr)
    errorFile.close()
    os.mkdir("output_"+px_nbr)
    patient_path = code_folder+"/output_"+px_nbr

    #we'll need the initial number of cells for the patient, and this value does not change depending on 
    #treatment and replicates
    #the idea is to count the number of lines in the CSV file associated to the patient 
    csvFileName = 'patientsCSV/'+px_IC
    rowCount = 0
    for row in open(csvFileName):
        rowCount += 1
    initialCells = rowCount

    #loop through your simulations
    #i is the simulation index
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

        #you will have to go get the file and copy it into "cells.csv " in the PhysiCell folder AND in the config folder
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
            newV0 = ((5/9740)*initialCells*initialCells)/(20*20*math.pi)

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

        """
        #Compress the files into a .zip
        os.chdir(patient_path)
        zip_output=px_nbr+'_'+file_name_main+'.zip'
        with zipfile.ZipFile(zip_output, 'w') as f:
            for file in glob.glob(treatment_folder+'/*'):
                f.write(file)

        now = datetime.now()
        date_string = now.strftime("%d/%m/%Y %H:%M")
        errorFile = open(errorFilePath,"a")
        errorFile.write("\n"+date_string+": Starting the transfer of the file "+zip_output)
        errorFile.write('\n')
        errorFile.close()

        #Send the file via Wetransfer
        #wt = WeTransfer(sender="blanche.mongeon@umontreal.ca",
        #                receivers=["blanche.mongeon@umontreal.ca"],
        #                channel='',
        #                message='File for '+px_nbr+' - Treatment: '+file_name_main,
        #                progress=True,
        #)
        wt = WeTransfer()
        upload = wt.start(zip_output)
        wt.uploadFile(zip_output)
        errorFile = open(errorFilePath,"a")
        errorFile.write("\n"+upload)
        errorFile.write('\n')

        now = datetime.now()
        date_string = now.strftime("%d/%m/%Y %H:%M")
        errorFile = open(errorFilePath,"a")
        errorFile.write("\n"+date_string+": Transfer of the file "+zip_output+" completed.")
        errorFile.write('\n')
        errorFile.close()
        

        #Delete the output file
        #os.chdir(code_folder)
        #shutil.rmtree(treatment_folder,ignore_errors=True,onerror=None)
        """
    errorFile = open(errorFilePath,"a")
    errorFile.write("\n "+px_nbr+" completed!")
    errorFile.close()

errorFile = open(errorFilePath,"a")
errorFile.write("\n All patients completed!")
errorFile.close()