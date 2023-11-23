# Treating glioblastoma treatment with immunotherapies: Spatial computational modelling illuminates the role of the tumour microenvironment
Author of the README file: Blanche Mongeon, Université de Montréal, blanche.mongeon@umontreal.ca

## Basic information about the code and how to use it

The **patientsCSV** folder contains the CSV files with the cell specifications for the bottomRepresentative (no adaptive immune system) and the topRepresentative. 
The **simulationsRepresentatives.py** file in the **python_scripts** folder executes 3 runs of each treatment combination analysed (none, TMZ, ICI, OV, OV+ICI, OV+ICI+TMZ) for each of the virtual patients. To run the file, you can use the command line
```
python3.10 simulations.py project
```
with **python_scripts** as your working directory. You have to make sure that python/3.10.1 is loaded before.

### If you want to run the code with a new CSV file
Little needs to be modified in the simulationsRepresentatives.py file in that case. Just make sure that the **patientsCSV** folder contains the CSV file with your cell specifications. Each row is a cell and needs to have 4 columns: 

column 1: x-position

column 2: y-position

column 3: z-position **this is set to 0 as we are simulating in 2D**

column 4: cell type (1=CD4 T cell, 2=Tumor cell, 3=CD8 T cell, 4=Stroma, 5=Macrophage)

If this is saved as **patientX.csv**, then on **line 62** of the simulationsRepresentatives.py file, add **patientX** to the array (do not include the ".csv"). If you don't want to run the simulations for the two virtual representatives, you can delete them from the array. 

### Using MATLAB to visualise results
The **matlab** folder contains a few function that might be useful to analyse the output files of a simulation. Refer to http://www.mathcancer.org/blog/working-with-physicell-snapshots-in-matlab/ to learn how to use them. 