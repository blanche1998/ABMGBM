function [cells,substrates,infectedCancerCells,internalizedVirus] = processMATLAB(patient,tx_folder,nOutput,nReplicates)
    
    %nOutput = str2num(nOutput);
    %nReplicates = str2num(nReplicates);
    tx = split(tx_folder,'/');
    tx = tx{end};
    cell_matrix_name = strcat(patient,'_',tx,'_cells');
    substrates_matrix_name = strcat(patient,'_',tx,'_substrates');
    rep_template = '/rep';
    %output_folder = strcat(tx,'/',rep_template,int2str(1));
    %MCDS = read_MultiCellDS_xml('initial.xml',output_folder);
    x_len = 200;%size(MCDS.mesh.X_coordinates,2);
    y_len = 200;%size(MCDS.mesh.Y_coordinates,2);
    nC = 5; %number of cell types
    %nS = MCDS.continuum_variables(:);
    nS = 6;%size(nS,1); %number of substrates
    cells = zeros(nC,nReplicates,nOutput+1);
    substrates = zeros(nS,nReplicates,nOutput+1,x_len,y_len);
    infectedCancerCells = zeros(nReplicates,nOutput+1);
    internalizedVirus = zeros(nReplicates,nOutput+1);

    
    num_output = '0000000';

    %loop through the replicates
    for n_rep=1:nReplicates
        output_folder = strcat(tx_folder,rep_template,int2str(n_rep));
        for k=0:nOutput
            num_output_temp = strcat(num_output,int2str(k));
            num_output_temp = num_output_temp(end-7:end);
    
            filename = strcat('output',num_output_temp,'.xml');
            
            MCDStemp = read_MultiCellDS_xml(filename,output_folder);
            types = MCDStemp.discrete_cells.metadata.type;
            
            dead_indexes = MCDStemp.discrete_cells.dead_cells;

            %count the number of cells of each type 
            for nType=1:nC
                indType = find(types==nType);
                indType = setdiff(indType,dead_indexes); %remove dead cells of that type
                nbr = size(indType,2);
                cells(nType,n_rep,k+1) = nbr;
            end
            
            %fill the substrate matrix
            for nSub=1:nS
                substrates(nSub,n_rep,k+1,:,:) = MCDStemp.continuum_variables(nSub).data(:,:);
            end

            %infection dynamics
            indexes_cancer = find(types==2); %indexes of cancer cells
            indexes_cancer = setdiff(indexes_cancer,dead_indexes); %remove dead cancer cells
            virus_intra_amounts = MCDStemp.discrete_cells.custom.intracellular_virus_amount(1,indexes_cancer);
            nstar = 1;%10; %infection threshold
            indexes_infected = find(virus_intra_amounts>nstar);
            %keep only the internalized amounts of infected cells
            virus_intra_amounts = virus_intra_amounts(1,indexes_infected); 
            infectedCancerCells(n_rep,k+1) = size(indexes_infected,2);
            internalizedVirus(n_rep,k+1) = sum(virus_intra_amounts,'all');

        end
    end
    %save everything
    %save(cell_matrix_name,"cells")
    %save(substrates_matrix_name,"substrates")
    
end

%}