# Treating glioblastoma treatment with immunotherapies: Spatial computational modelling illuminates the role of the tumour microenvironment
Author of the README file: Blanche Mongeon, Université de Montréal, blanche.mongeon@umontreal.ca

## Description of data files

Time course data are over 7 days and contain 85 hourly timepoints, from 0 to 168 hours in increments of 2 hours. Data files are .mat and can be opened with Matlab. Patients' indexes are ordered as: B1, B2, B3, T1, T2, T3. Treatment indexes are ordered as No treatment, TMZ, ICI, OV, TMZ+ICI+OV, ICI+OV. Cell types are ordered as CD4+ T cell, Cancer, CD8+ T cell, Stroma cells and Macrophages. 

**cancerPatients.mat** contain cancer cell counts for all 6 patients (B1, B2, B3, T1, T2, and T3) over 7 days depending on the simulated treatment. It is a 6x6x10x85 array; cancerPatients(i,j,k,l) is the l<sup>th</sup> timepoint for the k<sup>th</sup> replicate of patient j for the i<sup>th</sup> treatment regimen. For example, cancerPatients(2,2,5,4) is the the cancer cell count after 8 hours for the 5th replicate of patient B2 when TMZ is administered in monotherapy. Similarly, **cancerVirtualPatients.mat** contain cancer cell counts for our two virtual patients (BottomRepresentative and TopRepresentative) over 7 days depending on the simulated treatment. It is a 6x2x10x85 array; cancerPatients(i,j,k,l) is the l<sup>th</sup> timepoint for the k<sup>th</sup> replicate of virtual patient j for the i<sup>th</sup> treatment regimen.

**initialCellDistributionsPatients.mat** contain the initial cell distributions for all 6 patients. It is a 6x5 array, where initialCellDistributionsPatients(i,j) contains the fraction of initial cells for patient i that are of type j. Similarly, **initialCellDistributionsVirtualPatients.mat** contain the initial cell distributions for our two virtual patients. It is a 2x5 array, where initialCellDistributionsPatients(i,j) contains the fraction of initial cells for patient i that are of type j. 

**characteristics.mat** contains the data from panel A in Figure 3 _not min-max normalized_. To get the min-max normalized value, one can follow the following pseudocode:
```
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
```

The correlations between the characteristics (Figure 3 Panel B) can be obtained from the characteristics table once it is min-max normalized.

**cancerSA.mat** contains the results from the sensitivity analysis performed to analysis the effect of the cell types on treatment response (Panel A, figure 4). The results are _not_ min-max normalized; cancerSA(i,j) contains the fold change in cancer cells after 7 days when (j-1)*1000 cells of type i are added to the initial configuration of 5000 cells, with a TMZ+ICI+OV simulated treatment. 

**virionsSA.mat** contains the time course of virions when we increase the cancer cell count during our sensitivity analysis (Panel B, figure 4). virionsSA is a 6x85x3 array. virionsSA(i,j,k) is the virus density after 2*(j-1) hours of treatment initialisation (TMZ+ICI+OV) for the k<sup>th</sup> replicate when 1000*(i-1) supplementary cancer cells were added to the initial 5000 total cell configuration. 
