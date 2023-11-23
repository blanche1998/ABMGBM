clearvars
close all 

%% Define paths
%Change this to the folder with the output directories 
%simulationFolder = '/Volumes/Richard/Recherche/npjSBAPerspective/simulations3_6Px';
simulationFolder = '/NOBACKUP/mongeonb/PerspectiveReviewGBM/PhysiCell/Code/simulationsResults';
%patients_simul.csv file should be in the matlab folder
matlab_folder = cd(simulationFolder);
cd(matlab_folder)
visualisation_folder = strcat(simulationFolder,'/figures');

%% Define metrics

%load patients
%patients = load('patients_simul.mat');
%patients = patients.patients_simul;

%for i=1:length(patients)
%    temp = patients{i};
%    temp = convertStringsToChars(temp);
%    temp = temp(1:end-4); %remove '.csv' to keep only the patients' names
%    patients{i} = temp;
%end

%patients = {'CellsCSV'};
patients = {'bottomRepCell','topRepCell'};
Npx = 2; %number of patients

tx = {'None','TMZ','ICI','OV','TMZ+ICI+OV','ICI+OV'};%
Nsimulation = 6;
cellsTypes = ["TH","Cancer","CTL","Stroma","Macrophages"];
cont_variables = {'oxygen','TMZ','wall','chemokine','ICI','virus'};

Noutputs = 252;%36; %number of output files per patient
time = [0:Noutputs].*2;%in hours

Nreplicate = 3;

type_TH = 1;
type_cancer =2;
type_CTL = 3;
type_stroma = 4;
type_macrophage = 5;
NcellsTypes = 5;

cellsAll = zeros(Nsimulation,NcellsTypes,Npx,Nreplicate,length(time));
InfectedCells = zeros(Npx,Nsimulation,Nreplicate,length(time));
VirusAmounts = zeros(Npx,Nsimulation,Nreplicate,length(time));
xlen=200;
ylen=200;
%substratesAll = zeros(Nsimulation,length(cont_variables),Npx,Nreplicate,length(time),xlen,ylen);
%oxygenAll = zeros(Nsimulation,Npx,Nreplicate,length(time),xlen,ylen);
%TMZAll = zeros(Nsimulation,Npx,Nreplicate,length(time),xlen,ylen);
%wallAll = zeros(Nsimulation,Npx,Nreplicate,length(time),xlen,ylen);
%chemokineAll = zeros(Nsimulation,Npx,Nreplicate,length(time),xlen,ylen);
%ICIAll = zeros(Nsimulation,Npx,Nreplicate,length(time),xlen,ylen);
%virusAll = zeros(Nsimulation,Npx,Nreplicate,length(time),xlen,ylen);

folds = zeros(Nsimulation,Npx,NcellsTypes);

%initial immune proportion
immune_prop_initial = zeros(Npx,1);
%final immune proportion
immune_prop_final = zeros(Npx,Nsimulation);

special = zeros(0,3); %will keep track of special cases of increase 
                      % (when initial=0, but final!=0)

current = zeros(1,2);


%load('simulations6.mat')

%% loop through patients
for n=Npx:Npx
    current(1,1) = n;
    %save('current.mat','current')
    px = patients{n};
    tx1 = strcat(simulationFolder,'/output_',px,'/None');
    tx2 = strcat(simulationFolder,'/output_',px,'/TMZ');
    tx3 = strcat(simulationFolder,'/output_',px,'/ICI');
    tx4= strcat(simulationFolder,'/output_',px,'/OV');
    tx5 = strcat(simulationFolder,'/output_',px,'/TMZ+ICI+OV');
    tx6 = strcat(simulationFolder,'/output_',px,'/ICI+OV');

    
    folders = {tx1,tx2,tx3,tx4,tx5,tx6};%
    
    %we should already be in that file, but we want to make sure
    cd(matlab_folder)
    
    for nFolder=1:Nsimulation
        current(1,2) = nFolder;
        save('current.mat','current')%keep track of where you are in the simulations

        [cells_temp,~,infected_temp,virus_temp] = processMATLAB(convertCharsToStrings(px),...
            convertCharsToStrings(folders{nFolder}),Noutputs,Nreplicate);
        cellsAll(nFolder,:,n,:,:) = cells_temp;
        InfectedCells(n,nFolder,:,:) = infected_temp;
        VirusAmounts(n,nFolder,:,:) = virus_temp;
        %substratesAll(nFolder,:,n,:,:,:) = substrates_temp;
        %oxygenAll(nFolder,n,:,:,:,:) = substrates_temp(1,:,:,:,:);
        %TMZAll(nFolder,n,:,:,:,:) = substrates_temp(2,:,:,:,:);
        %wallAll(nFolder,n,:,:,:,:) = substrates_temp(3,:,:,:,:);
        %chemokineAll(nFolder,n,:,:,:,:) = substrates_temp(4,:,:,:,:);
        %ICIAll(nFolder,n,:,:,:,:) = substrates_temp(5,:,:,:,:);
        %virusAll(nFolder,n,:,:,:,:) = substrates_temp(6,:,:,:,:);

        initial_cells = mean(squeeze(cells_temp(:,:,1)),2)';
        immune_prop_initial(n,1) = (initial_cells(1,type_TH)+initial_cells(1,type_CTL))./sum(initial_cells);
        final_cells = mean(squeeze(cells_temp(:,:,end)),2)';
        immune_prop_final(n,nFolder) = (final_cells(1,type_TH)+final_cells(1,type_CTL))./sum(final_cells);
        %we have special cases and need to loop through each of the cells
        for nn=1:length(initial_cells)
            %nn will be the cell type index
            initial = initial_cells(1,nn);
            final = final_cells(1,nn);
            if(final==initial)
                folds(nFolder,n,nn) = 0;
            elseif(initial==0)
                %in that case final!=0 and thus we have an increase but
                %can't calculate a fold increase
                folds(nFolder,n,nn) = final;
                special(end+1,:) = [nFolder,n,nn];
            else
                inc = ((final/initial)-1)*100;%in percent
                folds(nFolder,n,nn) = inc;
            end
        end
        save('simulationsRepresentatives21days20230831.mat')
    end

end

%don't forget to save the workspace afterwards
save('simulationsRepresentatives21days20230831.mat')

