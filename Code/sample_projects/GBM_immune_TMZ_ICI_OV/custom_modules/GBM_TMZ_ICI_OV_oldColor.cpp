/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./GBM_TMZ_ICI_OV.h"
#include "../modules/PhysiCell_settings.h"

#include <filesystem>
#include <cmath>
#include <iostream>
#include <random>
#include <fstream>

Cell_Definition* cancer_cell;
Cell_Definition* Macrophage;
Cell_Definition* TH_cell;
Cell_Definition* stroma_cell;
Cell_Definition* CTL_cell;

void create_TH_cells( void )
{
	TH_cell = find_cell_definition( "TH cell" );
		
	// proliferation 
	TH_cell->functions.cycle_model = Ki67_basic;
	TH_cell->phenotype.cycle.sync_to_cycle_model( Ki67_basic); 
	int cycle_start_index = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_negative ); 
	int cycle_end_index = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_positive ); 
	TH_cell->phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = parameters.doubles("TH_prolif_rate"); 	
	TH_cell->phenotype.cycle.data.transition_rate(cycle_end_index,cycle_start_index) = parameters.doubles("TH_quiescent_transistion_rate"); 

	// cell actions
	TH_cell->phenotype.secretion.uptake_rates[1] = 0.0;
	TH_cell->phenotype.secretion.secretion_rates[1] = 0.0;
	TH_cell->phenotype.motility.migration_speed = parameters.ints("TH_migration_speed");
	
	// cell morphology
	TH_cell->phenotype.geometry.radius = 3.6;
	TH_cell->phenotype.volume.total = 185.66;
	TH_cell->phenotype.volume.fluid_fraction = 0.75;
	TH_cell->phenotype.volume.fluid = TH_cell->phenotype.volume.fluid_fraction*TH_cell->phenotype.volume.total;
	TH_cell->phenotype.volume.solid = TH_cell->phenotype.volume.total-TH_cell->phenotype.volume.fluid;
	TH_cell->phenotype.volume.nuclear = 95.21;
	TH_cell->phenotype.volume.nuclear_solid = 23.8;
	TH_cell->phenotype.volume.nuclear_fluid = TH_cell->phenotype.volume.nuclear - TH_cell->phenotype.volume.nuclear_solid;
	TH_cell->phenotype.volume.cytoplasmic = TH_cell->phenotype.volume.total - TH_cell->phenotype.volume.nuclear;
	TH_cell->phenotype.volume.cytoplasmic_fluid = TH_cell->phenotype.volume.fluid_fraction*TH_cell->phenotype.volume.cytoplasmic;
	TH_cell->phenotype.volume.cytoplasmic_solid = TH_cell->phenotype.volume.cytoplasmic-TH_cell->phenotype.volume.cytoplasmic_fluid;
	TH_cell->phenotype.volume.cytoplasmic_to_nuclear_ratio = 1.05;
	TH_cell->phenotype.volume.target_solid_cytoplasmic = TH_cell->phenotype.volume.cytoplasmic_solid;
	TH_cell->phenotype.volume.target_solid_nuclear = TH_cell->phenotype.volume.nuclear_solid;
	TH_cell->phenotype.volume.target_fluid_fraction = TH_cell->phenotype.volume.fluid_fraction;
	TH_cell->phenotype.volume.target_cytoplasmic_to_nuclear_ratio = TH_cell->phenotype.volume.cytoplasmic_to_nuclear_ratio;
	
	// update function
	TH_cell->functions.update_phenotype = TH_functions;
	
	return;
}

void create_CTL_cells( void )
{
	CTL_cell = find_cell_definition( "CTL cell" );
	
	// proliferation 
	CTL_cell->functions.cycle_model = Ki67_basic;
	CTL_cell->phenotype.cycle.sync_to_cycle_model( Ki67_basic); 
	int cycle_start_index = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_negative ); 
	int cycle_end_index = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_positive); 
	CTL_cell->phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = parameters.doubles("CTL_prolif_rate"); 
	CTL_cell->phenotype.cycle.data.transition_rate(cycle_end_index,cycle_start_index) = parameters.doubles("CTL_quiescent_transistion_rate");	
	
	// cell actions
	CTL_cell->phenotype.secretion.uptake_rates[1] = 0.0;
	CTL_cell->phenotype.secretion.secretion_rates[1] = 0.0;
	CTL_cell->phenotype.motility.migration_speed = parameters.ints("CTL_migration_speed");
	
	//cell morphology
	CTL_cell->phenotype.geometry.radius = 3.6;
	CTL_cell->phenotype.volume.total = 185.66;
	CTL_cell->phenotype.volume.fluid_fraction = 0.75;
	CTL_cell->phenotype.volume.fluid = CTL_cell->phenotype.volume.fluid_fraction*CTL_cell->phenotype.volume.total;
	CTL_cell->phenotype.volume.solid = CTL_cell->phenotype.volume.total-CTL_cell->phenotype.volume.fluid;
	CTL_cell->phenotype.volume.nuclear = 96.23;
	CTL_cell->phenotype.volume.nuclear_solid = 24.06;
	CTL_cell->phenotype.volume.nuclear_fluid = CTL_cell->phenotype.volume.nuclear - CTL_cell->phenotype.volume.nuclear_solid;
	CTL_cell->phenotype.volume.cytoplasmic = CTL_cell->phenotype.volume.total - CTL_cell->phenotype.volume.nuclear;
	CTL_cell->phenotype.volume.cytoplasmic_fluid = CTL_cell->phenotype.volume.fluid_fraction*CTL_cell->phenotype.volume.cytoplasmic;
	CTL_cell->phenotype.volume.cytoplasmic_solid = CTL_cell->phenotype.volume.cytoplasmic-CTL_cell->phenotype.volume.cytoplasmic_fluid;
	CTL_cell->phenotype.volume.cytoplasmic_to_nuclear_ratio = 1.03;
	CTL_cell->phenotype.volume.target_solid_cytoplasmic = CTL_cell->phenotype.volume.cytoplasmic_solid;
	CTL_cell->phenotype.volume.target_solid_nuclear = CTL_cell->phenotype.volume.nuclear_solid;
	CTL_cell->phenotype.volume.target_fluid_fraction = CTL_cell->phenotype.volume.fluid_fraction;
	CTL_cell->phenotype.volume.target_cytoplasmic_to_nuclear_ratio = CTL_cell->phenotype.volume.cytoplasmic_to_nuclear_ratio;
	
	// cell update phenotype		
	CTL_cell->functions.update_phenotype = CTL_functions;
	
	return;
}

void create_stroma_cells( void )
{
	stroma_cell = find_cell_definition( "stroma cell" );
	
	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live); 
	stroma_cell->phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = 0; 

	// cell actions
	static int TMZ_index = microenvironment.find_density_index( "TMZ");
	stroma_cell->phenotype.secretion.uptake_rates[TMZ_index] = parameters.doubles("stroma_virus_uptake_rate")*0.5;
	stroma_cell->phenotype.secretion.secretion_rates[1] = 0.0;
	stroma_cell->phenotype.motility.migration_speed = 0;
	stroma_cell->phenotype.motility.is_motile = false;
	stroma_cell->phenotype.molecular.fraction_released_at_death[TMZ_index] = 0;
	
	// cell morphology
	stroma_cell->phenotype.geometry.radius = parameters.doubles("stroma_radius");
	stroma_cell->phenotype.volume.total = 1767;
	stroma_cell->phenotype.volume.fluid_fraction = 0.75;
	stroma_cell->phenotype.volume.fluid = stroma_cell->phenotype.volume.fluid_fraction*stroma_cell->phenotype.volume.total;
	stroma_cell->phenotype.volume.solid = stroma_cell->phenotype.volume.total-stroma_cell->phenotype.volume.fluid;
	stroma_cell->phenotype.volume.nuclear = 500;
	stroma_cell->phenotype.volume.nuclear_solid = 125;
	stroma_cell->phenotype.volume.nuclear_fluid = stroma_cell->phenotype.volume.nuclear - stroma_cell->phenotype.volume.nuclear_solid;
	stroma_cell->phenotype.volume.cytoplasmic = stroma_cell->phenotype.volume.total - stroma_cell->phenotype.volume.nuclear;
	stroma_cell->phenotype.volume.cytoplasmic_fluid = stroma_cell->phenotype.volume.fluid_fraction*stroma_cell->phenotype.volume.cytoplasmic;
	stroma_cell->phenotype.volume.cytoplasmic_solid = stroma_cell->phenotype.volume.cytoplasmic-stroma_cell->phenotype.volume.cytoplasmic_fluid;
	stroma_cell->phenotype.volume.cytoplasmic_to_nuclear_ratio = 2.53;
	stroma_cell->phenotype.volume.target_solid_cytoplasmic = stroma_cell->phenotype.volume.cytoplasmic_solid;
	stroma_cell->phenotype.volume.target_solid_nuclear = stroma_cell->phenotype.volume.nuclear_solid;
	stroma_cell->phenotype.volume.target_fluid_fraction = stroma_cell->phenotype.volume.fluid_fraction;
	stroma_cell->phenotype.volume.target_cytoplasmic_to_nuclear_ratio = stroma_cell->phenotype.volume.cytoplasmic_to_nuclear_ratio;
	
	// update phenotype
	stroma_cell->functions.update_phenotype = stroma_function;
	
	return;
}

void create_macrophages( void )
{
	Macrophage = find_cell_definition( "Macrophage" ); 
		
	// proliferation 
	Macrophage->functions.cycle_model = Ki67_basic;
	Macrophage->phenotype.cycle.sync_to_cycle_model( Ki67_basic); 
	int cycle_start_index = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_negative ); 
	int cycle_end_index = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_positive); 
	Macrophage->phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = parameters.doubles("CTL_prolif_rate"); 
	Macrophage->phenotype.cycle.data.transition_rate(cycle_end_index,cycle_start_index) = parameters.doubles("CTL_quiescent_transistion_rate"); 	
	
	// cell actions
	Macrophage->phenotype.secretion.uptake_rates[1] = 0.0;
	Macrophage->phenotype.secretion.secretion_rates[1] = 0.0;
	Macrophage->phenotype.motility.migration_speed = parameters.ints("CTL_migration_speed");
	
	//cell morphology
	Macrophage->phenotype.geometry.radius = 10.5;
	Macrophage->phenotype.volume.total = 4849;
	Macrophage->phenotype.volume.fluid_fraction = 0.75;
	Macrophage->phenotype.volume.fluid = Macrophage->phenotype.volume.fluid_fraction*Macrophage->phenotype.volume.total;
	Macrophage->phenotype.volume.solid = Macrophage->phenotype.volume.total-Macrophage->phenotype.volume.fluid;
	Macrophage->phenotype.volume.nuclear = 485;
	Macrophage->phenotype.volume.nuclear_solid = 121.25;
	Macrophage->phenotype.volume.nuclear_fluid = Macrophage->phenotype.volume.nuclear - Macrophage->phenotype.volume.nuclear_solid;
	Macrophage->phenotype.volume.cytoplasmic = Macrophage->phenotype.volume.total - Macrophage->phenotype.volume.nuclear;
	Macrophage->phenotype.volume.cytoplasmic_fluid = Macrophage->phenotype.volume.fluid_fraction*Macrophage->phenotype.volume.cytoplasmic;
	Macrophage->phenotype.volume.cytoplasmic_solid = Macrophage->phenotype.volume.cytoplasmic-Macrophage->phenotype.volume.cytoplasmic_fluid;
	Macrophage->phenotype.volume.cytoplasmic_to_nuclear_ratio = 9;
	Macrophage->phenotype.volume.target_solid_cytoplasmic = Macrophage->phenotype.volume.cytoplasmic_solid;
	Macrophage->phenotype.volume.target_solid_nuclear = Macrophage->phenotype.volume.nuclear_solid;
	Macrophage->phenotype.volume.target_fluid_fraction = Macrophage->phenotype.volume.fluid_fraction;
	Macrophage->phenotype.volume.target_cytoplasmic_to_nuclear_ratio = Macrophage->phenotype.volume.cytoplasmic_to_nuclear_ratio;
	
	// cell update phenotype		
	Macrophage->functions.update_phenotype = Macrophage_functions;
	
	return;
}


void create_cancer_cells( void )
{
	cancer_cell = find_cell_definition( "cancer cell" ); 
	
	// cell actions - TMZ
	static int TMZ_index = microenvironment.find_density_index( "TMZ");
	cancer_cell->phenotype.secretion.uptake_rates[TMZ_index] = parameters.doubles("cell_TMZ_uptake_rate"); 
	cancer_cell->phenotype.secretion.secretion_rates[TMZ_index] = 0.0;
	cancer_cell->phenotype.motility.migration_speed = 0.0;
	// cell actions - OV (taken from Adrianne's code)
	static int virus_index = microenvironment.find_density_index( "virus");
	cancer_cell->phenotype.secretion.uptake_rates[virus_index] = 0.0;
	cancer_cell->phenotype.secretion.secretion_rates[virus_index] = 0.0;
	cancer_cell->phenotype.motility.migration_speed = 0.0;
	cancer_cell->phenotype.molecular.fraction_released_at_death[virus_index] = 1;
	
	// cell morphology
	cancer_cell->phenotype.geometry.radius = parameters.doubles("R_cell_GBM");
	cancer_cell->phenotype.volume.total = 4/3*3.1416*(cancer_cell->phenotype.geometry.radius)*(cancer_cell->phenotype.geometry.radius)*(cancer_cell->phenotype.geometry.radius);//5203.7;
	cancer_cell->phenotype.volume.fluid_fraction = parameters.doubles("f_F");
	cancer_cell->phenotype.volume.fluid = cancer_cell->phenotype.volume.fluid_fraction*cancer_cell->phenotype.volume.total;
	cancer_cell->phenotype.volume.solid = cancer_cell->phenotype.volume.total-cancer_cell->phenotype.volume.fluid;
	cancer_cell->phenotype.volume.nuclear = parameters.doubles("V_N_GBM");;
	cancer_cell->phenotype.volume.nuclear_solid = 185;
	cancer_cell->phenotype.volume.nuclear_fluid = cancer_cell->phenotype.volume.nuclear - cancer_cell->phenotype.volume.nuclear_solid;
	cancer_cell->phenotype.volume.cytoplasmic = cancer_cell->phenotype.volume.total - cancer_cell->phenotype.volume.nuclear;
	cancer_cell->phenotype.volume.cytoplasmic_fluid = cancer_cell->phenotype.volume.fluid_fraction*cancer_cell->phenotype.volume.cytoplasmic;
	cancer_cell->phenotype.volume.cytoplasmic_solid = cancer_cell->phenotype.volume.cytoplasmic-cancer_cell->phenotype.volume.cytoplasmic_fluid;
	cancer_cell->phenotype.volume.cytoplasmic_to_nuclear_ratio = cancer_cell->phenotype.volume.cytoplasmic/cancer_cell->phenotype.volume.nuclear;
	cancer_cell->phenotype.volume.target_solid_cytoplasmic = cancer_cell->phenotype.volume.cytoplasmic_solid;
	cancer_cell->phenotype.volume.target_solid_nuclear = cancer_cell->phenotype.volume.nuclear_solid;
	cancer_cell->phenotype.volume.target_fluid_fraction = cancer_cell->phenotype.volume.fluid_fraction;
	cancer_cell->phenotype.volume.target_cytoplasmic_to_nuclear_ratio = cancer_cell->phenotype.volume.cytoplasmic_to_nuclear_ratio;
	cancer_cell->phenotype.volume.calcified_fraction = 0; 
	cancer_cell->phenotype.volume.calcification_rate = 0;	
	
	// set functions 
	cancer_cell->functions.update_phenotype = cancer_cell_proliferation_infection_movement; 
	
	return; 
}

void create_cell_types( void )
{
	// set the random seed 
	SeedRandom( parameters.ints("random_seed") );  
		
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );
	cell_defaults.phenotype.molecular.sync_to_microenvironment( &microenvironment );
	
	// Make sure we're ready for 2D
	cell_defaults.functions.set_orientation = up_orientation;  
	cell_defaults.phenotype.geometry.polarity = 1.0; 
	cell_defaults.phenotype.motility.restrict_to_2D = true; 
	
	//setting cycle model to live
	cell_defaults.phenotype.cycle.sync_to_cycle_model( live ); 
	
	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
	
	cell_defaults.phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) = 0;
	
 	// turn off death
	int apoptosis_index = cell_defaults.phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model ); 
	cell_defaults.phenotype.death.rates[apoptosis_index] = 0.0; 

	// reduce cell velocity
	cell_defaults.phenotype.motility.migration_speed = 0.0;
	cell_defaults.custom_data.add_variable( "persistence_time", "dimensionless", 0.0 ); 
	cell_defaults.custom_data.add_variable( "cell_motility_type", "dimensionless", 0.0 );
	cell_defaults.custom_data.add_variable( "rep_rate", "dimensionless", 0.0 );
	cell_defaults.custom_data.add_variable( "attachment lifetime" , "min" , 0 ); 
	cell_defaults.custom_data.add_variable( "Tumour_antigen_expression","dimensionless",0.0); 
	cell_defaults.custom_data.add_variable( "Activation boolean","dimensionless",0.0);
	cell_defaults.custom_data.add_variable( "Attachment boolean","dimensionless",0.0);
	cell_defaults.custom_data.add_variable( "time_to_next_phagocytosis",0.0);
	
	// add variables to track virus infection start time, length of time and amount of virus
	cell_defaults.custom_data.add_variable( "intracellular_virus_amount", "dimensionless", 0.0 ); // amount of intracellular virus
	//cell_defaults.custom_data.add_variable( "persistence_time", "dimensionless", 0.0 ); // how long cells will persist in move or stop phenotype
	//cell_defaults.custom_data.add_variable( "cell_motility_type", "dimensionless", 0.0 );
	//cell_defaults.custom_data.add_variable( "attachment lifetime" , "min" , 0 ); // how long it can stay attached 

	Parameter<double> paramD;
	
	paramD = parameters.doubles[ "elastic_coefficient" ]; 
	cell_defaults.custom_data.add_variable( "elastic coefficient" , paramD.units, paramD.value ); 

	// turn off secretion from these cells (oxygen)
	cell_defaults.phenotype.secretion.secretion_rates[0] = 0;
	cell_defaults.phenotype.secretion.uptake_rates[0] = 0;
	cell_defaults.phenotype.secretion.saturation_densities[0] = 10; 

	static int virus_index = microenvironment.find_density_index( "virus");
	static int wall_index = microenvironment.find_density_index( "wall");
	static int chemokine_index = microenvironment.find_density_index( "chemokine");
	static int TMZ_index = microenvironment.find_density_index( "TMZ");
	
	cell_defaults.phenotype.secretion.secretion_rates[virus_index] = 0;
	cell_defaults.phenotype.secretion.uptake_rates[virus_index] = 0.0;
	cell_defaults.phenotype.secretion.saturation_densities[virus_index] = parameters.doubles("rhostar_virus");//;//10; 
	cell_defaults.phenotype.molecular.fraction_released_at_death[virus_index] = 0;//1;

	cell_defaults.phenotype.molecular.fraction_released_at_death[wall_index] = 1;//1;
	cell_defaults.phenotype.secretion.secretion_rates[wall_index] = 0;
	cell_defaults.phenotype.secretion.uptake_rates[wall_index] = 0;
	cell_defaults.phenotype.secretion.saturation_densities[wall_index] = 10; 
	
	cell_defaults.phenotype.molecular.fraction_released_at_death[chemokine_index] = 1;//1;
	cell_defaults.phenotype.secretion.secretion_rates[chemokine_index] = 0;
	cell_defaults.phenotype.secretion.uptake_rates[chemokine_index] = 0;
	cell_defaults.phenotype.secretion.saturation_densities[chemokine_index] = parameters.doubles("chemokine_saturation_density");//5; 
	
	cell_defaults.phenotype.secretion.secretion_rates[TMZ_index] = 0;
	cell_defaults.phenotype.secretion.uptake_rates[TMZ_index] = 0.0;
	cell_defaults.phenotype.secretion.saturation_densities[TMZ_index] = parameters.doubles("virus_saturation_density");
	cell_defaults.phenotype.molecular.fraction_released_at_death[TMZ_index] = 1;
	cell_defaults.phenotype.molecular.fraction_released_at_death[0] = 1;
	cell_defaults.phenotype.molecular.fraction_released_at_death[2] = 1;
	cell_defaults.phenotype.molecular.fraction_released_at_death[3] = 1;
	
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 
	
	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	cell_defaults.functions.update_phenotype = cancer_cell_proliferation_infection_movement;
	cell_defaults.phenotype.motility.is_motile = true; 
	
	//HERE_BM
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
	
	// create the cell types
	create_cancer_cells();
	create_macrophages();
	create_TH_cells();
	create_stroma_cells();
	create_CTL_cells();
		
	build_cell_definitions_maps(); 
	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	
	// set domain parameters 
	
	// make sure not override and go back to 2D 
	if( default_microenvironment_options.simulate_2D == false )
	{
		std::cout << "Warning: overriding XML config option and setting to 2D!" << std::endl; 
		default_microenvironment_options.simulate_2D = true; 
	}
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment();
	
	static int TMZ_index = microenvironment.find_density_index( "TMZ");
	static int oxygen_index = microenvironment.find_density_index( "oxygen");
	static int wall_index = microenvironment.find_density_index( "wall");
	static int chemokine_index = microenvironment.find_density_index( "chemokine");
	static int ICI_index = microenvironment.find_density_index( "ICI");
	static int virus_index = microenvironment.find_density_index( "virus");
	
	for( int n = 0 ; n < microenvironment.mesh.voxels.size(); n++ )
	{	
		std::vector<double> ECMdense = microenvironment.mesh.voxels[n].center; 
		//assign random ECM density to microenvironment voxel
		
		if( ECMdense[0]*ECMdense[0]+ECMdense[1]*ECMdense[1]>(parameters.doubles("tumour_radius")+10)*(parameters.doubles("tumour_radius")+10))
		{	
			microenvironment(n)[oxygen_index] = 0;
			microenvironment(n)[TMZ_index] = 0;
			microenvironment(n)[wall_index] = 1;
			microenvironment(n)[chemokine_index] = 0;
			microenvironment(n)[ICI_index] = 0;
			microenvironment(n)[virus_index] = 0;
			
		}
		else if(ECMdense[0]*ECMdense[0]+ECMdense[1]*ECMdense[1]>(parameters.doubles("tumour_radius")-10)*(parameters.doubles("tumour_radius")-10))
		{
			
			microenvironment(n)[oxygen_index] = 0;
			microenvironment(n)[TMZ_index] = 0;
			microenvironment(n)[wall_index] = 3.5;
			microenvironment(n)[chemokine_index] = 0;
			microenvironment(n)[ICI_index] = 0;
			microenvironment(n)[virus_index] = 0;
		}
		else if( ECMdense[0]*ECMdense[0]+ECMdense[1]*ECMdense[1]>250*250 )
		{	
			microenvironment(n)[oxygen_index] = 0;
			microenvironment(n)[TMZ_index] = 0;
			microenvironment(n)[wall_index] = 5;
			microenvironment(n)[chemokine_index] = 0;
			microenvironment(n)[ICI_index] = 0;
			microenvironment(n)[virus_index] = 0;
		}
		else
		{   
			microenvironment(n)[oxygen_index] = 0;
			microenvironment(n)[TMZ_index] = 0;
			microenvironment(n)[wall_index] = 10;
			microenvironment(n)[chemokine_index] = 0;
			microenvironment(n)[ICI_index] = 0;
			microenvironment(n)[virus_index] = 0;
		}
			
	}

	//the OV_tx is a boolean parameter that allows us to change the administered treatment
	if ( parameters.bools("OV_tx") ) {
		//this could be added as a parameter to the Physicell setting file 
		std::string injection = "unique";
		adminOV(injection,virus_index);
	} else {
		std::cout << "we don't give any OV.";
	}

	return; 
}

//administration of the initial injection(s) of the OV 
void adminOV( std::string typeInjection , int ind )
{
	if (typeInjection == "unique") {

		

		for( int n = 0 ; n < microenvironment.mesh.voxels.size(); n++ )
		{
			std::vector<double> ECMdense = microenvironment.mesh.voxels[n].center;
			if( (ECMdense[0])*(ECMdense[0])+(ECMdense[1])*(ECMdense[1])<20*20)//10*10)// Centre
			{
				microenvironment(n)[ind] = 7.12;//parameters.doubles("V0");//7.12; //PFU*C0/(pi*20^2) 
			}
		}
		return;
	} else if (typeInjection == "multiple") {

		//this is for one administration but multiple injections
		// locations of dense patches
		double x_offset_patch_1 = 173.38;
		double y_offset_patch_1 = -270.5;
		double x_offset_patch_2 = 41.32;
		double y_offset_patch_2 = 1015.1;
		double x_offset_patch_3 = -757.84;
		double y_offset_patch_3 = 128.4;
		double x_offset_patch_4 = 899.79;
		double y_offset_patch_4 = -391.6;
		double x_offset_patch_5 = -9.99;
		double y_offset_patch_5 = -859.9;

		for( int n = 0 ; n < microenvironment.mesh.voxels.size(); n++ )
		{
			std::vector<double> ECMdense = microenvironment.mesh.voxels[n].center; 
			if( (x_offset_patch_1-ECMdense[0])*(x_offset_patch_1-ECMdense[0])+(y_offset_patch_1-ECMdense[1])*(y_offset_patch_1-ECMdense[1])<20*20)// Centre of patch 1
			{
				microenvironment(n)[ind] = 7.12;
			}
			else if((x_offset_patch_2-ECMdense[0])*(x_offset_patch_2-ECMdense[0])+(y_offset_patch_2-ECMdense[1])*(y_offset_patch_2-ECMdense[1])<20*20)// Centre of patch 2
			{
				microenvironment(n)[ind] = 7.12;
			}
			else if((x_offset_patch_3-ECMdense[0])*(x_offset_patch_3-ECMdense[0])+(y_offset_patch_3-ECMdense[1])*(y_offset_patch_3-ECMdense[1])<20*20) //Centre of patch 3
			{
				microenvironment(n)[ind] = 7.12;
			}
			else if((x_offset_patch_4-ECMdense[0])*(x_offset_patch_4-ECMdense[0])+(y_offset_patch_4-ECMdense[1])*(y_offset_patch_4-ECMdense[1])<20*20) // Centre of patch 4
			{	
				microenvironment(n)[ind] = 7.12;
			}
			else if((x_offset_patch_5-ECMdense[0])*(x_offset_patch_5-ECMdense[0])+(y_offset_patch_5-ECMdense[1])*(y_offset_patch_5-ECMdense[1])<20*20) // Centre of patch 5
			{
				microenvironment(n)[ind] = 7.12;
			}
			
		}
	} else {
		std::cout << "In the .cpp custom file, line 457, the OV injection can either be \"multiple\" or \"unique\". \n";
		std::cout << "Please enter the type of OV injections you want to administer: \n ";
		std::string newType;
		std::cin >> newType;
		adminOV(newType,ind);
	}
	return;
}

// function to make sure that we are injecting into tissue
std::vector<double> inject_tissue() {

	/*
	pugi::xml_document doc;
	pugi::xml_parse_result result = doc.load_file("../config/PhysiCell_settings_TMZ_ICI_OV.xml");
	pugi::xml_node root = doc.document_element();
	pugi::xml_node node = root.child( "initial_conditions" ); 
	node = root.child( "cell_positions" ); 

	// get filename 

	std::string folder = xml_get_string_value( node, "folder" ); 
	std::string filename = xml_get_string_value( node, "filename" ); 
	// std::string input_filename = folder + "/" + filename; 
	std::string filetype = node.attribute("type").value() ; 
	*/
	std::string input_filename = "cells.csv";

	//std::ifstream file( input_filename, std::ios::in );
	std::ifstream file;
	file.open( input_filename,std::ios::in );
	
	if( !file )
	{ 
		std::cout << "Error: " << input_filename << " not found during cell loading. Quitting." << std::endl; 
		exit(-1);
	}
	
	std::string line;
	double mini = 10000;
	std::vector<double> center = {0,0};
	while (std::getline(file, line))
	{
		//std::cout << "\n we're reading the file \n";
		std::vector<double> data;
		csv_to_vector( line.c_str() , data ); 

		std::vector<double> position = { data[0] , data[1] , data[2] };

		if ((position[0]*position[0]+position[1]*position[1]) < 20*20) {
			file.close();
			center[0] = position[0];
			center[1] = position[1];
			//std::filesystem::current_path("path"); //setting path
			return center;
		} else {
			double x_current = position[0];
			double y_current = position[1];
			if (x_current*x_current + y_current*y_current < mini) {
				mini = x_current*x_current + y_current*y_current;
				center[0] = position[0];
				center[1] = position[1];
			} 
		}

	}

	file.close(); 
	//std::filesystem::current_path("path"); //setting path
	return center;
}

//administration of the initial injection(s) of the OV with a different center
void adminOV( std::string typeInjection , int ind, double x_center, double y_center )
{
	if ( parameters.bools("OV_tx") ) {
		std::cout << "\n we'll be injecting the virus with center " << x_center << "," << y_center << "\n";
	
		for( int n = 0 ; n < microenvironment.mesh.voxels.size(); n++ )
		{
			std::vector<double> ECMdense = microenvironment.mesh.voxels[n].center;
			double r_x = ECMdense[0] - x_center;
			double r_y = ECMdense[1] - y_center;
			if ( (ECMdense[0])*(ECMdense[0])+(ECMdense[1])*(ECMdense[1])<10*10) {
				std::cout << "\n" << n;
				microenvironment(n)[ind] = 0; //we reset it to 0
			}
			//problem here
			if( ((r_x)*(r_x)+(r_y)*(r_y))<20*20)// Centre
			{
				//std::cout << n;
				microenvironment(n)[ind] = 7.12;//parameters.doubles("V0");//7.12;
			}
		}

	}
	
	return;
}

void setup_tissue_circle_immune( void )
{
	double Radius = parameters.doubles("tumour_radius");
	
	double x = 0.0;
	double y = 0.0;
	
	// setting up distributions for movement and persistance of cells
	std::vector<double> go_times_cumul(8);
    go_times_cumul[0] = 0.01;
	go_times_cumul[1] = 0.962;
	go_times_cumul[2] = 0.9735;
	go_times_cumul[3] = 0.9835;
	go_times_cumul[4] = 0.9935;
	go_times_cumul[5] = 0.9955;
	go_times_cumul[6] = 0.9975;
	go_times_cumul[7] = 1;
	
	std::vector<double> persistence_times_vec(8);
    persistence_times_vec[0] = 0;
	persistence_times_vec[1] = 30;
	persistence_times_vec[2] = 60;
	persistence_times_vec[3] = 90;
	persistence_times_vec[4] = 120;
	persistence_times_vec[5] = 150;
	persistence_times_vec[6] = 180;
	persistence_times_vec[7] = 240;
	
	std::vector<double> speed_cumul(12);
    speed_cumul[0] = 0.0014;
	speed_cumul[1] = 0.0317;
	speed_cumul[2] = 0.2441;
	speed_cumul[3] = 0.5137;
	speed_cumul[4] = 0.7598;
	speed_cumul[5] = 0.8822;
	speed_cumul[6] = 0.9453;
	speed_cumul[7] = 0.9787;
	speed_cumul[8] = 0.9882;
	speed_cumul[9] = 0.9937;
	speed_cumul[10] = 0.9963;
	speed_cumul[11] = 1;
	
	std::vector<double> speed_vec(12);
    speed_vec[0] = 0.0833;
	speed_vec[1] = 0.1667;
	speed_vec[2] = 0.25;
	speed_vec[3] = 0.333;
	speed_vec[4] = 0.4167;
	speed_vec[5] = 0.5;
	speed_vec[6] = 0.5833;
	speed_vec[7] = 0.667;
	speed_vec[8] = 0.75;
	speed_vec[9] = 0.833;
	speed_vec[10] = 0.9167;
	speed_vec[11] = 1;
		
		
	double GBM_NO = parameters.ints("initial_GBM_cells");
	double stroma_NO = parameters.ints("initial_stroma_cells");
	std::default_random_engine generator;

	double nu_mean = parameters.doubles("virus_replication_rate");
	double nu_variance = 0.01;
	double shape = nu_mean*nu_mean/nu_variance;
	double scale = nu_variance/nu_mean;
	std::gamma_distribution<double> distribution(shape,scale);	
	
	load_cells_from_pugixml(); 	
	
	for( int n=0; n < (*all_cells).size(); n++ )
	{
		Cell* pC = (*all_cells)[n]; 
		if(pC->type ==2)
		{
		int persistence_time_index = pC->custom_data.find_variable_index( "persistence_time" );
		int cell_motility_type_index = pC->custom_data.find_variable_index( "cell_motility_type" );
		
		static int Tumour_antigen_expression_index = pC->custom_data.find_variable_index( "Tumour_antigen_expression" );
		pC->custom_data.variables[Tumour_antigen_expression_index].value = 1;
		
		double p = UniformRandom();
		if(p<=0.5)// GO
		{
			pC->custom_data.variables[cell_motility_type_index].value = 1;
			double speed_var = UniformRandom();
			
			for( int k=0; k<12; )
			{
				if( speed_var> speed_cumul[k] )
				{k++;}
				else
				{
					pC->phenotype.motility.migration_speed = speed_vec[k];
					k = 12;
				}
			}
		} 
		else
		{pC->custom_data.variables[cell_motility_type_index].value = 2;} // STOP
	
		double go_stop_var = UniformRandom();
		
		for( int j=0; j<8; )
		{
			if( go_stop_var> go_times_cumul[j] )
			{j++;}
			else
			{
				pC->custom_data.variables[persistence_time_index].value = persistence_times_vec[j];
				j = 8;
			}
			
		}	
		}
	}
	
	return;
	
}

void cancer_cell_proliferation_infection_movement( Cell* pCell, Phenotype& phenotype, double dt )
{

	double R = parameters.doubles("R_cell_GBM");
	double SA = 4*3.1416*R*R;
	double xi_val = parameters.doubles("maximum_cell_density");
	double s =sqrt(2.0/(sqrt(3.0)*xi_val));
	double pressure = 6*(1-((1/(2*R))*s))*(1-((1/(2*R))*s));
	double pressure_scale = 0.027288820670331; 
	double max_pressure = pressure/pressure_scale;	

	//tumour cell proliferation
	if(pCell->type ==2)
	{		
		static int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
		static int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
		
		static int apoptosis_index = pCell->phenotype.death.find_death_model_index( "apoptosis" );
		
		double EC50 = parameters.doubles("IC50");
	
		static int TMZ_index = microenvironment.find_density_index("TMZ");	
		double TMZ = pCell->nearest_density_vector()[TMZ_index];
		double TMZ_effect = TMZ/(TMZ+EC50);
			
		pCell->phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) = parameters.doubles("GBM_cell_proliferation_rate")*(1-TMZ_effect);		
		
		pCell->phenotype.death.rates[apoptosis_index] = parameters.doubles("GBM_cell_proliferation_rate")*(TMZ_effect);
				
		if( pCell->state.simple_pressure>max_pressure) // if cell under too much pressure -> no proliferation
		{pCell->phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) = 0;}
		

	}

	
	
	cell_movement( pCell, phenotype, dt);
	
	//this can happen if we're injecting an oncolytic virus
	if ( parameters.bools("OV_tx") ) {
		infection_dynamics( pCell, phenotype, dt );
	}
	return;
	
}

void cell_movement( Cell* pCell, Phenotype& phenotype, double dt ) 
{
	static int wall_index = microenvironment.find_density_index( "wall" ); 
	double wall_amount = pCell->nearest_density_vector()[wall_index];
	
	static int persistence_time_index = pCell->custom_data.find_variable_index( "persistence_time" );
	static int cell_motility_type_index = pCell->custom_data.find_variable_index( "cell_motility_type" );
	
	double persistence_time = pCell->custom_data.variables[persistence_time_index].value;
	double cell_motility_type = pCell->custom_data.variables[cell_motility_type_index].value; 
	
	std::vector<double> go_times_cumul(8);
    go_times_cumul[0] = 0.01;
	go_times_cumul[1] = 0.962;
	go_times_cumul[2] = 0.9735;
	go_times_cumul[3] = 0.9835;
	go_times_cumul[4] = 0.9935;
	go_times_cumul[5] = 0.9955;
	go_times_cumul[6] = 0.9975;
	go_times_cumul[7] = 1;
	
	std::vector<double> persistence_times_vec(8);
    persistence_times_vec[0] = 0;
	persistence_times_vec[1] = 30;
	persistence_times_vec[2] = 60;
	persistence_times_vec[3] = 90;
	persistence_times_vec[4] = 120;
	persistence_times_vec[5] = 150;
	persistence_times_vec[6] = 180;
	persistence_times_vec[7] = 240;
	
	std::vector<double> speed_cumul(12);
    speed_cumul[0] = 0.0014;
	speed_cumul[1] = 0.0317;
	speed_cumul[2] = 0.2441;
	speed_cumul[3] = 0.5137;
	speed_cumul[4] = 0.7598;
	speed_cumul[5] = 0.8822;
	speed_cumul[6] = 0.9453;
	speed_cumul[7] = 0.9787;
	speed_cumul[8] = 0.9882;
	speed_cumul[9] = 0.9937;
	speed_cumul[10] = 0.9963;
	speed_cumul[11] = 1;
	
	std::vector<double> speed_vec(12);
    speed_vec[0] = 0.0833;
	speed_vec[1] = 0.1667;
	speed_vec[2] = 0.25;
	speed_vec[3] = 0.333;
	speed_vec[4] = 0.4167;
	speed_vec[5] = 0.5;
	speed_vec[6] = 0.5833;
	speed_vec[7] = 0.667;
	speed_vec[8] = 0.75;
	speed_vec[9] = 0.833;
	speed_vec[10] = 0.9167;
	speed_vec[11] = 1;
	
	if( wall_amount<2 & pCell->type == 2 )
	{
		pCell->phenotype.motility.migration_speed = 0;
	}
	else if(pCell->type != 2)
	{
		pCell->phenotype.motility.migration_speed = 4;
	}
	else		
	{
		if( persistence_time <= PhysiCell_globals.current_time ) // if the cell's persistence time is up
		{
			// assign new type (stop = 2, or go = 1)
			double new_type_rand = UniformRandom();
			if(new_type_rand<=0.5)// GO
			{
				pCell->custom_data.variables[cell_motility_type_index].value = 1; // assign go type
				
				double speed_var = UniformRandom();
			
				for( int k=0; k<12; )
				{
					if( speed_var> speed_cumul[k] )
					{k++;}
					else
					{
						pCell->phenotype.motility.migration_speed = speed_vec[k]; // assign migration speed
						k = 12;
					}
				}
			} 
			else
			{pCell->custom_data.variables[cell_motility_type_index].value = 2;
			pCell->phenotype.motility.migration_speed = 0;} // assign STOP type

			// assign persistence time - needs to be a real time!
			double go_stop_var = UniformRandom();
			for( int j=0; j<8; )
			{
				if( go_stop_var> go_times_cumul[j] )
				{j++;}
					else
				{
					pCell->custom_data.variables[persistence_time_index].value = persistence_times_vec[j]+PhysiCell_globals.current_time; // assign persist time
					j = 8;
				}
			}
		}
			
	}
	return;
}

void infection_dynamics( Cell* pCell, Phenotype& phenotype, double dt )
{
			
	//cell infection
	static int virus_signal_index = microenvironment.find_density_index( "virus");
	static int apoptosis_model_index =	pCell->phenotype.death.find_death_model_index( "apoptosis" );
	static double intracellular_virus_index = pCell->custom_data.find_variable_index( "intracellular_virus_amount" );
		
	double n = pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index];//custom_data.variables[intracellular_virus_index].value;
	double p = pCell->nearest_density_vector()[virus_signal_index];
	double u = parameters.doubles("u_g");
	double Vvoxel = microenvironment.mesh.voxels[0].volume;//volume of voxel
	double nstar = parameters.doubles("m_half");//10;//infection threshol
	double alp = parameters.doubles("alpha");//1000;//virus burst number
	double nu = parameters.doubles("gamma");
	double pmax = parameters.doubles("rho_max");//0.0125;
	
	if( n>1)
	{
	//	std::cout<<"m_i = "<<n<<" rho_virus = "<<p<<" u_g = "<<u<<" mstar = "<<nstar<<" alph = "<<alp<<" gamma = "<<nu<<" rho_max = "<<pmax<<std::endl;
	}
	pCell->custom_data.variables[intracellular_virus_index].value = pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index];		
		
	if( pCell->phenotype.death.dead == false )// cell not dead	
	{	
		if(p<pmax)
		{pCell->phenotype.secretion.uptake_rates[virus_signal_index] = u*p/(n/Vvoxel+nstar/Vvoxel);}
		else
		{pCell->phenotype.secretion.uptake_rates[virus_signal_index] = u*pmax*pmax/(n/Vvoxel+nstar/Vvoxel)/p;}
		
		if( n > 1 && n <= alp) // update amount inside due to replication 
		{
			pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index] = n+dt*(nu*n); 
			pCell->custom_data.variables[intracellular_virus_index].value = pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index];
		}
		else if( n <= 1 )
		{
			pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index] = 0;
			pCell->custom_data.variables[intracellular_virus_index].value = pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index];
		}
		else if( n > alp-1)
		{
			pCell->custom_data.variables[intracellular_virus_index].value = pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index];
			//pCell->phenotype.molecular.fraction_released_at_death[virus_signal_index] = 1;
			if( pCell->type==2 )
			{
				pCell->phenotype.molecular.fraction_released_at_death[virus_signal_index] = 1;
			}
			pCell->phenotype.secretion.uptake_rates[virus_signal_index] = 0;
			pCell->start_death( apoptosis_model_index );	
		}
		else if( pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index]<0 )
		{std::cout<<"Negative intracellular virus!! m_i = "<<pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index]<<std::endl;}
	}
	if( pCell->phenotype.death.dead == true && pCell->phenotype.molecular.fraction_released_at_death[virus_signal_index]>0)
	{
		virus_induced_lysis(pCell, phenotype, dt );
	}
	
	return;
}

void virus_induced_lysis( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int virus_signal_index = microenvironment.find_density_index( "virus");
	double pstar = pCell->phenotype.secretion.saturation_densities[virus_signal_index];
	double delta_V = parameters.doubles("delta_V");//0.1466;
	double Vvoxel = microenvironment.mesh.voxels[0].volume;//volume of voxel
	double p = pCell->nearest_density_vector()[virus_signal_index];
	double n = pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index];
	
	if( pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index]> 1 ) 
	{
			double amount_to_add = (n-n*exp(-delta_V*dt))/Vvoxel;
			if( amount_to_add > pstar-p )
			{
				pCell->nearest_density_vector()[virus_signal_index] += (pstar-p)*dt;
				pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index] -= (pstar-p)*Vvoxel*dt;
			}
			else
			{
				pCell->nearest_density_vector()[virus_signal_index] += (amount_to_add)*dt;
				pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index] = n*exp(-delta_V*dt);
			}	
	}

	return;
}

void TMZ_induced_death( Cell* pCell, Phenotype& phenotype, double dt ) 
{
	
	static int TMZ_index = microenvironment.find_density_index( "TMZ");
	static int apoptosis_model_index =	pCell->phenotype.death.find_death_model_index( "apoptosis" );
	static int Activation_boolean_index = pCell->custom_data.find_variable_index("Activation boolean");
	
	double TMZ = pCell->nearest_density_vector()[TMZ_index];
	double IC50 = parameters.doubles("IC50");
	
	double Imax = 1*(TMZ/(TMZ+IC50));
	double pd = Imax;//1-exp(Imax*dt);
		
	if( UniformRandom()<pd)
	{
		pCell->start_death( apoptosis_model_index );
		pCell->phenotype.molecular.fraction_released_at_death[TMZ_index] = 0;
		if(pCell->state.neighbors.size()>0)
		{dettach_cells( pCell, pCell->state.neighbors[0] );} 
	}
	return;
	
}

void dettach_cells( Cell* pCell_1 , Cell* pCell_2 )
{
	#pragma omp critical
	{
		bool found = false; 
		int i = 0; 
		while( !found && i < pCell_1->state.neighbors.size() )
		{
			// if cell 2 is in cell 1's list, remove it
			if( pCell_1->state.neighbors[i] == pCell_2 )
			{
				int n = pCell_1->state.neighbors.size(); 
				// copy last entry to current position 
				pCell_1->state.neighbors[i] = pCell_1->state.neighbors[n-1]; 
				// shrink by one 
				pCell_1->state.neighbors.pop_back(); 
				found = true; 
			}
			i++; 
		}
	
		found = false; 
		i = 0; 
		while( !found && i < pCell_2->state.neighbors.size() )
		{
			// if cell 1 is in cell 2's list, remove it
			if( pCell_2->state.neighbors[i] == pCell_1 )
			{
				int n = pCell_2->state.neighbors.size(); 
				// copy last entry to current position 
				pCell_2->state.neighbors[i] = pCell_2->state.neighbors[n-1]; 
				// shrink by one 
				pCell_2->state.neighbors.pop_back(); 
				found = true; 
			}
			i++; 
		}

	}
	
	return; 
}

void Macrophage_functions( Cell* pCell, Phenotype& phenotype, double dt )
{

	static int apoptosis_index = phenotype.death.find_death_model_index( "Apoptosis" ); 
	static int Activation_boolean_index = pCell->custom_data.find_variable_index("Activation boolean");

	if( pCell->custom_data[Activation_boolean_index] < 1 )
	{ pCell->phenotype.death.rates[apoptosis_index] = 0; }
	else
	{ pCell->phenotype.death.rates[apoptosis_index] = parameters.doubles("normal_mac_death_rate");}

	if( pCell->phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL;
		pCell->functions.custom_cell_rule = NULL; 
		return; 
	}

	// check for cells to eat 
	std::vector<Cell*> neighbors = pCell->cells_in_my_container(); 

	if( neighbors.size() < 2 )
	{ return; } 
	static int CD8_Tcell_type = 3; 
	static int CD4_Tcell_type = 1;  
	static int cancer_type = 2;  
	int n = 0; 
	Cell* pContactCell = neighbors[n]; 
	static int virus_signal_index = microenvironment.find_density_index( "virus");
	double nstar = parameters.doubles("infection_threshold");
	while( n < neighbors.size() )
	{
		pContactCell = neighbors[n]; 
		
		double cell_cell_distance = sqrt((pContactCell->position[0]-pCell->position[0])*(pContactCell->position[0]-pCell->position[0])+(pContactCell->position[1]-pCell->position[1])*(pContactCell->position[1]-pCell->position[1]));
		double radius_mac = pCell->phenotype.geometry.radius; 
		double radius_test_cell = pContactCell->phenotype.geometry.radius; 
		if( pContactCell != pCell && pContactCell->phenotype.death.dead == false && pContactCell->type == CD8_Tcell_type 
			&& pCell->custom_data["Activation boolean"] > 0 && cell_cell_distance<=parameters.doubles("epsilon_distance")*(radius_mac+radius_test_cell)) 
		{
			pContactCell->custom_data["Activation boolean"]=1;
			n=neighbors.size();
		}
		else if( pContactCell != pCell && pContactCell->phenotype.death.dead == false && pContactCell->type == CD4_Tcell_type 
			&& pCell->custom_data["Activation boolean"]> 0 && cell_cell_distance<=parameters.doubles("epsilon_distance")*(radius_mac+radius_test_cell))
		{
			
			pContactCell->custom_data["Activation boolean"]=1;
			n=neighbors.size();
		} 	
		n++;
	}
	if( pCell->phenotype.volume.total> parameters.doubles("threshold_macrophage_volume"))
	{
		pCell->phenotype.death.rates[apoptosis_index] = parameters.doubles("exhausted_macrophage_death_rate");
		return;
	}	
	int time_to_next_phagocytosis_index = pCell->custom_data.find_variable_index( "time_to_next_phagocytosis" );
	if( pCell->custom_data.variables[time_to_next_phagocytosis_index].value>PhysiCell_globals.current_time )
	{return;}	
	
	double probability_of_phagocytosis = 1;
	double material_internalisation_rate = parameters.doubles("material_internalisation_rate"); 

	static int chemokine_index = microenvironment.find_density_index("chemokine" ); 
	// static int virus_signal_index = microenvironment.find_density_index( "virus");
	// double nstar = parameters.doubles("infection_threshold");
	
	n = 0; 
	Cell* pTestCell = neighbors[n]; 
	while( n < neighbors.size() )
	{
		pTestCell = neighbors[n]; 
		if( pTestCell != pCell && pTestCell->phenotype.death.dead == true &&  
			UniformRandom() < probability_of_phagocytosis ) 
		{
			
			{
				double volume_ingested_cell = pTestCell->phenotype.volume.total;
				pCell->ingest_cell( pTestCell ); 
				double time_to_ingest = volume_ingested_cell*material_internalisation_rate;
				pCell->custom_data.variables[time_to_next_phagocytosis_index].value = PhysiCell_globals.current_time+time_to_ingest;				
			}	
			pCell->phenotype.motility.migration_speed = parameters.doubles("activated_speed"); 
			pCell->custom_data.variables[Activation_boolean_index].value = 1.0; 
			pCell->phenotype.secretion.secretion_rates[chemokine_index] = parameters.doubles("chemokine_secretion_rate");	
			
			return; 
		} else if (pTestCell != pCell && pTestCell->phenotype.death.dead == false && pTestCell->phenotype.molecular.internalized_total_substrates[virus_signal_index] > nstar ) {
				// if the cell being tested is not myself and the cell is infected with some virus above a given threshold
	
				pCell->custom_data.variables[Activation_boolean_index].value = 1.0; 
				// macrophage becomes activated
				return; 
		} 	
		n++; 
	}			
	return; 
}


void macrophage_mechanics( Cell* pCell, Phenotype& phenotype, double dt )
{

	if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL;
		pCell->functions.custom_cell_rule = NULL; 

		return; 
	}

	return; 
}

void TH_functions( Cell* pCell, Phenotype& phenotype, double dt )
{
	
	static int wall_index = microenvironment.find_density_index( "wall" ); 
	double wall_amount = pCell->nearest_density_vector()[wall_index];
	
	std::vector<double> ae_ini(3);	
	
	if( wall_amount<2 )
	{
		pCell->phenotype.motility.migration_speed = 4;
		ae_ini = -1*pCell->position;
		pCell->phenotype.motility.migration_bias = 1;
		normalize( &( ae_ini ) );
		pCell->phenotype.motility.migration_bias_direction = ae_ini;
	}
	else 
	{
		pCell->phenotype.motility.migration_speed = 4;
		pCell->phenotype.motility.migration_bias = 0;
	}
	std::vector<Cell*> nearby = pCell->cells_in_my_container();
	
	static int TMZ_index = microenvironment.find_density_index( "TMZ");
	static int chemokine_index = microenvironment.find_density_index( "chemokine");
	static int virus_signal_index = microenvironment.find_density_index( "virus");
	
	static int Activation_boolean_index = pCell->custom_data.find_variable_index( "Activation boolean" );
	
	int cycle_start_index = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_negative ); 
	int cycle_end_index = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_positive ); 
	
	Cell* pC = NULL;
	bool stop = false;
	int i=0;

	double nstar = parameters.doubles("infection_threshold");
	
	if(pCell->custom_data[Activation_boolean_index]<1)
	{
		while( !stop && i < nearby.size() )
		{
			pC = nearby[i];
			int Tumour_antigen_expression_index = pC->custom_data.find_variable_index( "Tumour_antigen_expression" ); 
		
			if( pC->custom_data[Tumour_antigen_expression_index]>0 &&	pC->phenotype.death.dead == false 
				&& pC != pCell && pC->type == 2 && UniformRandom()>0.5)
			{stop = true;}

			//the next 4 lines are taken from Adrianne's code
			if( pC->phenotype.molecular.internalized_total_substrates[virus_signal_index] > nstar &&
			pC->phenotype.death.dead == false &&
			pC != pCell && pC->type != 4)
			{ stop = true; }
			
			i++;
		
			if( stop == false )
				{ pC = NULL; }
		}
		
		if( pC )
		{
			pCell->phenotype.secretion.secretion_rates[chemokine_index] = parameters.doubles("chemokine_secretion_rate");	
			pCell->phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = parameters.doubles("TH_prolif_rate")*parameters.doubles("TH_prolif_increase_due_to_stimulus"); 	
			pCell->phenotype.motility.migration_speed = 0.1;
			pCell->custom_data.variables[Activation_boolean_index].value=1;
		
		}
		else
		{
			pCell->phenotype.secretion.secretion_rates[chemokine_index] = 0;//0.1;	
				
			pCell->phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = parameters.doubles("TH_prolif_rate"); 
		}
	}
	return;
}

void stroma_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	
	if (parameters.bools("OV_tx")) {
		static int virus_signal_index = microenvironment.find_density_index("virus");
		static double intracellular_virus_index = pCell->custom_data.find_variable_index( "intracellular_virus_amount" );
		pCell->custom_data.variables[intracellular_virus_index].value = pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index];
		if(pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index]>1e5)
		{
			pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index] = 0;
			pCell->custom_data.variables[intracellular_virus_index].value = 0;
		}
		else
		{
			pCell->custom_data.variables[intracellular_virus_index].value = pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index];
		}
	}
	
	return;
}

void CTL_functions( Cell* pCell, Phenotype& phenotype, double dt )
{	
	
	if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL; 
		return; 
	}	
	
	int cycle_start_index = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_negative ); 
	int cycle_end_index = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_positive ); 
	
	static int Activation_boolean_index = pCell->custom_data.find_variable_index( "Activation boolean" ); 
	static int attach_lifetime_i = pCell->custom_data.find_variable_index( "attachment lifetime" ); 
	static int Attachment_boolean_index = pCell->custom_data.find_variable_index( "Attachment boolean" ); 
	std::vector<Cell*> neighbors = pCell->cells_in_my_container(); 
	static int bound_PD1_ICI=pCell->custom_data.find_variable_index( "Bound_PD1_receptors_on_cell" );

	if( neighbors.size()>0 && pCell->custom_data[Activation_boolean_index]<1)
	{
		std::vector<Cell*> nearby = pCell->cells_in_my_container(); 
		int i = 0; 
		while( i < nearby.size() )
		{
			static int Activation_boolean_index_check = nearby[i]->custom_data.find_variable_index( "Activation boolean" ); 
			if( nearby[i]->custom_data[Activation_boolean_index_check]>0 )
			{
				pCell->custom_data.variables[Activation_boolean_index].value = 1;
				i = nearby.size();
			}
			i++; 
		}
	}
	
	if( pCell->state.neighbors.size()>0 && pCell->custom_data[Attachment_boolean_index] > 0 )
	{
		extra_elastic_attachment_mechanics( pCell, phenotype, dt );
		bool dettach_me = false; 
		if( immune_cell_attempt_apoptosis( pCell, pCell->state.neighbors[0], dt ) )
		{
			immune_cell_trigger_apoptosis( pCell, pCell->state.neighbors[0] ); 
			if(pCell->type==3)
			{
			static int nR_EB = pCell->custom_data.find_variable_index( "Bound_PD1_receptors_on_cell" );
			static int nV_external = microenvironment.find_density_index( "ICI" );	
			}
			dettach_me = true; 
		}
		 
		if( dettach_me )
		{
			
			dettach_cells( pCell, pCell->state.neighbors[0] ); 
			phenotype.motility.is_motile = true; 
			CTL_movement( pCell, phenotype, dt);
			pCell->phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = parameters.doubles("CTL_prolif_rate")*parameters.doubles("CTL_prolif_increase_due_to_stimulus"); 	
			pCell->custom_data.variables[Attachment_boolean_index].value = 0;
		}

		return; 
	}
	else if( pCell->custom_data[Activation_boolean_index] > 0 && pCell->custom_data[Attachment_boolean_index] <1 && pCell->custom_data[bound_PD1_ICI] > parameters.ints("PD1_ICI_threshold"))
	{
		
		if( immune_cell_check_neighbors_for_attachment( pCell , dt) )
		{			
			phenotype.motility.is_motile = false; 
		
			return; 
		}
	}
	
	phenotype.motility.is_motile = true; 
	CTL_movement( pCell, phenotype, dt);
		
	return;
}

void extra_elastic_attachment_mechanics( Cell* pCell, Phenotype& phenotype, double dt )
{
	for( int i=0; i < pCell->state.neighbors.size() ; i++ )
	{
		add_elastic_velocity( pCell, pCell->state.neighbors[i], pCell->custom_data["elastic coefficient"] ); 
	}

	return; 
}	

void add_elastic_velocity( Cell* pActingOn, Cell* pAttachedTo , double elastic_constant )
{
	std::vector<double> displacement = pAttachedTo->position - pActingOn->position; 
	axpy( &(pActingOn->velocity) , elastic_constant , displacement ); 
	
	return; 
}

Cell* immune_cell_check_neighbors_for_attachment( Cell* pAttacker , double dt )
{
	std::vector<Cell*> nearby = pAttacker->cells_in_my_container(); 
	int i = 0; 
	while( i < nearby.size() )
	{
		if( nearby[i] != pAttacker )
		{
			if( immune_cell_attempt_attachment( pAttacker, nearby[i] , dt ) )
			{ return nearby[i]; }
		}
		i++; 
	}
	
	return NULL; 
}

bool immune_cell_attempt_attachment( Cell* pAttacker, Cell* pTarget , double dt )
{
	double max_attachment_distance = parameters.doubles("max_attachment_distance"); 
	
	static int attach_lifetime_i = pAttacker->custom_data.find_variable_index( "attachment lifetime" ); 
	static int Tumour_antigen_expression_index = pTarget->custom_data.find_variable_index( "Tumour_antigen_expression" ); 
	static int Attachment_boolean_index = pAttacker->custom_data.find_variable_index( "Attachment boolean" ); 
	
	static int virus_signal_index = microenvironment.find_density_index( "virus");
	double nstar = parameters.doubles("m_half"); // how long the cell needs to attach for  the infected cell to be killed
	double internal_virus = pTarget->phenotype.molecular.internalized_total_substrates[virus_signal_index];
	
	double kill_time = parameters.doubles("time_to_kill_cell");
	
	//the OR condition is taken from Adrianne's code
	if( (pTarget->custom_data[Tumour_antigen_expression_index] > 0 || internal_virus > nstar ) && pTarget->phenotype.death.dead == false && pTarget->type!=4)
	{
		std::vector<double> displacement = pTarget->position - pAttacker->position;
		double distance_scale = norm( displacement ); 
		if( distance_scale > max_attachment_distance )
		{ return false; } 
	
		attach_cells( pAttacker, pTarget ); 
		pAttacker->custom_data.variables[attach_lifetime_i].value = PhysiCell_globals.current_time + kill_time;
		pAttacker->custom_data.variables[Attachment_boolean_index].value = 1;
		return true; 
	}

	return false; 
}

void CTL_movement( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int wall_index = microenvironment.find_density_index( "wall" ); 
	double wall_amount = pCell->nearest_density_vector()[wall_index];
	
	static int chemokine_index = microenvironment.find_density_index( "chemokine" ); 
	double chemokine_amount = pCell->nearest_density_vector()[chemokine_index];
	
	int cycle_start_index = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_negative ); 
	int cycle_end_index = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_positive ); 
	
	
	std::vector<double> ae_ini(3);
	if( wall_amount<2 )
	{
		pCell->phenotype.motility.migration_speed = 1;
		pCell->phenotype.motility.migration_bias = 1;
		pCell->phenotype.motility.migration_bias_direction = pCell->nearest_gradient(wall_index);
		ae_ini = -1*pCell->position;
	
		pCell->phenotype.motility.migration_bias = 1;
		normalize( &( ae_ini ) );
		pCell->phenotype.motility.migration_bias_direction = ae_ini;
		
		return;
	}
	else if(chemokine_amount>1e-8)
	{
		pCell->phenotype.motility.migration_speed = parameters.doubles("CTL_min_speed")+(parameters.doubles("CTL_max_speed")-parameters.doubles("CTL_min_speed"))*(chemokine_amount/(parameters.doubles("Chemokine_EC50")+chemokine_amount));
		pCell->phenotype.motility.migration_bias = parameters.doubles("CTL_chemokine_migration_bias");
		pCell->phenotype.motility.migration_bias_direction = pCell->nearest_gradient(chemokine_index);
		pCell->phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = parameters.doubles("CTL_prolif_rate");
		
		return;
	}
	else if(chemokine_amount>1e-3)
	{
		pCell->phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = parameters.doubles("CTL_prolif_rate")*parameters.doubles("CTL_prolif_increase_due_to_stimulus"); 	
		pCell->phenotype.motility.migration_speed = parameters.doubles("CTL_min_speed")+(parameters.doubles("CTL_max_speed")-parameters.doubles("CTL_min_speed"))*(chemokine_amount/(parameters.doubles("Chemokine_EC50")+chemokine_amount));
		pCell->phenotype.motility.migration_bias = parameters.doubles("CTL_chemokine_migration_bias");
		pCell->phenotype.motility.migration_bias_direction = pCell->nearest_gradient(chemokine_index);
	}
	else
	{
		pCell->phenotype.motility.migration_bias = 0;
		pCell->phenotype.motility.migration_speed = 4;
		pCell->phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = 7.9026*1e-5;
		return;
	}
}

bool immune_cell_attempt_apoptosis( Cell* pAttacker, Cell* pTarget, double dt )
{
	static int attach_lifetime_i = pAttacker->custom_data.find_variable_index( "attachment lifetime" ); 
	if( pAttacker->custom_data[attach_lifetime_i] < PhysiCell_globals.current_time )
	{ 
		return true; 
	}
	return false; 
}

bool immune_cell_trigger_apoptosis( Cell* pAttacker, Cell* pTarget )
{
	static int apoptosis_model_index = pTarget->phenotype.death.find_death_model_index( "apoptosis" );	
	
	static int TMZ_index = microenvironment.find_density_index( "TMZ" ); 
	if( pTarget->phenotype.death.dead == true )
	{ return false; }

	pTarget->start_death( apoptosis_model_index );
	pTarget->phenotype.molecular.fraction_released_at_death[TMZ_index] = 0;	

	return true; 
}

void receptor_dynamics_model( Cell* pCell, Phenotype& phenotype, double dt )
{

	static int nV_external = microenvironment.find_density_index( "ICI" ); 
	
	static int nV_internal = pCell->custom_data.find_variable_index( "ICI" );  
	
	static int nR_EU = pCell->custom_data.find_variable_index( "Unbound_PD1_receptors_on_cell" );
	static int nR_EB = pCell->custom_data.find_variable_index( "Bound_PD1_receptors_on_cell" );
	
	static int nR_bind = pCell->custom_data.find_variable_index( "ICI_binding_rate" ); 	
	static int nR_unbind = pCell->custom_data.find_variable_index( "ICI_unbinding_rate" ); 
	if( phenotype.death.dead == true )
	{ return; } 
	if( pCell->type != 3 )
	{ return; } 
	double newly_bound = phenotype.molecular.internalized_total_substrates[nV_external];  
	double excess_binding = newly_bound - pCell->custom_data[nR_EU]; 
	if( excess_binding > 0.0 )
	{
		newly_bound = pCell->custom_data[nR_EU]; 
		static double one_virion_to_density = 1.0 / microenvironment.mesh.dV; 
		#pragma omp critical 
		{
			pCell->nearest_density_vector()[nV_external] += excess_binding * one_virion_to_density; 
		}
	}
	phenotype.molecular.internalized_total_substrates[nV_external] = 0.0; 
	
	pCell->custom_data[nR_EB] += newly_bound; 

	pCell->custom_data[nR_EU] -= newly_bound; 
	
	double dR_EU = dt*pCell->custom_data[nR_unbind]*pCell->custom_data[nR_EB];
	if( dR_EU > pCell->custom_data[nR_EB] )
	{ dR_EU = pCell->custom_data[nR_EB]; }
	pCell->custom_data[nR_EB] -= dR_EU; 
	pCell->custom_data[nR_EU] += dR_EU; 
	
	static double one_virion_to_density = 1.0 / microenvironment.mesh.dV;
	#pragma omp critical 
	{
		pCell->nearest_density_vector()[nV_external] += dR_EU*one_virion_to_density; 
	}
	
	phenotype.secretion.uptake_rates[nV_external] = 
		pCell->custom_data[nR_bind] * pCell->custom_data[nR_EU]; 
	
	return; 
}

 void receptor_dynamics_model( double dt )
{
	#pragma omp parallel for 
	for( int n=0; n < (*all_cells).size() ; n++ )
	{
		Cell* pC = (*all_cells)[n]; 
		receptor_dynamics_model( pC, pC->phenotype , dt ); 	
	}
	
	return; 
} 

std::vector<std::string> colouring( Cell* pCell )
{

	std::vector< std::string > output( 4, "darkgrey" ); 
	
	static int v_index = microenvironment.find_density_index( "virus");
	static double intracellular_virus_index = pCell->custom_data.find_variable_index( "intracellular_virus_amount" );
	
	double p_min = 1;
	double p_max = parameters.doubles("alpha");
	
	double n_I = pCell->phenotype.molecular.internalized_total_substrates[v_index];
	if(n_I>10000 && pCell->type ==2)
	{std::cout<<" High viral infection m_i = "<<n_I<<std::endl;}
	
	//dead
	if(pCell->phenotype.death.dead == true)
	{
		output[0] = "rgb(255, 255, 224)";
		output[1] = "rgb(255, 255, 224)";
		output[2] = "rgb(255, 228, 181)";
		output[3] = "rgb(255, 228, 181)";
		
		return output; 
	}	
	else if(pCell->type==1) //TH
	{
		int oncoprotein = (int) round( (1.0/(p_max-p_min)) * (pCell->phenotype.molecular.internalized_total_substrates[v_index]-p_min) * 255.0 ); 
		char szTempString [128]; // ceates a character array that can store 128
		sprintf( szTempString , "rgb(%u,%u,%u)", oncoprotein, oncoprotein, 255-oncoprotein ); // puts oncoprotein, oncoprotein and 255-oncoprotein in place of u u u
			
			
		if(pCell-> phenotype.cycle.data.current_phase_index==0) { 
					output[0] = "orange";
					output[1] = "orange";
					output[2] = "coral";
					output[3] = "coral";
					return output; 
		} else if(pCell-> phenotype.cycle.data.current_phase_index==1) {	
					output[0] = "darkred";
					output[1] = "darkred";
					output[2] = "firebrick";
					output[3] = "firebrick";
		}
	}
	else if( pCell->type == 2) //Cancer cell
	{	
		if( n_I>1) {
			double p_min = 1;
			int oncoprotein1 = (int) round((210-139)*(1.0/(p_max-p_min)) * (n_I-1)); 
			int oncoprotein2 = (int) round((180-69)*(1.0/(p_max-p_min)) * (n_I-1)); 
			int oncoprotein3 = (int) round((140-19)*(1.0/(p_max-p_min)) * (n_I-1)); 
			
			char szTempString [128]; // ceates a character array that can store 128
			sprintf( szTempString , "rgb(%u,%u,%u)", 210-oncoprotein1, 180-oncoprotein2, 140-oncoprotein3); // puts oncoprotein, oncoprotein and 255-oncoprotein in place of u u u
			
			output[0].assign( szTempString );
			output[1]="brown";
			output[2].assign( szTempString );
			output[3]="brown";
			
			return output;
		} else {
			output[0] = "rgb(104, 55, 99)";//"orchid";//"rgb(255,230,230)";
			output[1] = "rgb(104, 55, 99)";
			output[2] = "rgb(85, 50, 70)";//"plum";//"rgb(255,230,230)";
			output[3] = "rgb(85, 50, 70)";

			return output; 
		}
	}
	else if( pCell->type == 3) //CTL
	{
		if(pCell-> phenotype.cycle.data.current_phase_index==0) {   
			output[0] = "aquamarine";
			output[1] = "lightsteelblue";
			output[2] = "lightskyblue";
			output[3] = "aquamarine";
		} else if(pCell-> phenotype.cycle.data.current_phase_index==1) {	
	
			output[0] = "darkslateblue";
			output[1] = "darkblue";
			output[2] = "darkblue";
			output[3] = "aquamarine";
		} 					
	}
	else if(pCell->type == 4) //STROMA
	{ 
		output[0] = "rgb(234, 172, 199)";
		output[1] = "rgb(234, 172, 199)";
		output[2] = "rgb(243, 186, 211)";
		output[3] = "rgb(243, 186, 211)";		
		return output; 
	}
	else if(pCell->type == 5) //MACROPHAGE
	{
		output[0] = "green";
		output[1] = "green";
		output[2] = "green";
		output[3] = "green";		
		return output; 
	}
	return output;
}

