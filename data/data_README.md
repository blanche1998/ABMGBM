# Treating glioblastoma treatment with immunotherapies: Spatial computational modelling illuminates the role of the tumour microenvironment
Author of the README file: Blanche Mongeon, Université de Montréal, blanche.mongeon@umontreal.ca

## Description of data files

Time course data are over 7 days and contain 85 hourly timepoints, from 0 to 168 hours in increments of 2 hours. Data files are .mat and can be opened with Matlab. Patients' indexes are ordered as: B1, B2, B3, T1, T2, T3. Treatment indexes are ordered as No treatment, TMZ, ICI, OV, TMZ+ICI+OV, ICI+OV. Cell types are ordered as CD4+ T cell, Cancer, CD8+ T cell, Stroma cells and Macrophages. 

**cancerPatients.mat** contain cancer cell counts for all 6 patients (B1, B2, B3, T1, T2, and T3) over 7 days depending on the simulated treatment. It is a 6x6x10x85 array; cancerPatients(i,j,k,l) is the l<sup>th</sup> timepoint for the k<sup>th</sup> replicate of patient j for the i<sup>th</sup> treatment regimen. For example, cancerPatients(2,2,5,4) is the the cancer cell count after 8 hours for the 5th replicate of patient B2 when TMZ is administered in monotherapy. Similarly, **cancerVirtualPatients.mat** contain cancer cell counts for our two virtual patients (BottomRepresentative and TopRepresentative) over 7 days depending on the simulated treatment. It is a 6x2x10x85 array; cancerPatients(i,j,k,l) is the l<sup>th</sup> timepoint for the k<sup>th</sup> replicate of virtual patient j for the i<sup>th</sup> treatment regimen.

**initialCellDistributionsPatients.mat** contain the initial cell distributions for all 6 patients. It is a 6x5 array, where initialCellDistributionsPatients(i,j) contains the fraction of initial cells for patient i that are of type j. Similarly, **initialCellDistributionsVirtualPatients.mat** contain the initial cell distributions for our two virtual patients. It is a 2x5 array, where initialCellDistributionsPatients(i,j) contains the fraction of initial cells for patient i that are of type j. 

**characteristics.mat** contains the data from panel A in Figure 3 _not min-max normalized_. To get the min-max normalized value, 

maxR = size(characteristics,1);
maxC = size(characteristics,2);
for j=1:maxC
    max = max(characteristics(:,j)); 
    min = min(characteristics(:,j));
    for i=1:maxR
        currentVal = characteristics(i,j);
        normVal = (currentVal-min)./(max-min);
        characteristics(i,j) = normVal;
    end for
end for



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