/*
#############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the ver-  #
# sion number, such as below:                                               #
#                                                                           #
# We implemented and solved the model using PhysiCell (Version 1.0.0) [1].  #
#                                                                           #
# [1] A Ghaffarizadeh, SH Friedman, SM Mumenthaler, and P Macklin,          #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for            #
#     Multicellular Systems, 2016 (in preparation).                         #
#                                                                           #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite       #
#     BioFVM as below:                                                      #
#                                                                           #
# We implemented and solved the model using PhysiCell (Version 1.0.0) [1],  #
# with BioFVM [2] to solve the transport equations.                         #
#                                                                           #
# [1] A Ghaffarizadeh, SH Friedman, SM Mumenthaler, and P Macklin,          #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for            #
#     Multicellular Systems, 2016 (in preparation).                         #
#                                                                           #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient     #
#    parallelized diffusive transport solver for 3-D biological simulations,#
#    Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730 #
#                                                                           #
#############################################################################
#                                                                           #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)   #
#                                                                           #
# Copyright (c) 2015-2016, Paul Macklin and the PhysiCell Project           #
# All rights reserved.                                                      #
#                                                                           #
# Redistribution and use in source and binary forms, with or without        #
# modification, are permitted provided that the following conditions are    #
# met:                                                                      #
#                                                                           #
# 1. Redistributions of source code must retain the above copyright notice, #
# this list of conditions and the following disclaimer.                     #
#                                                                           #
# 2. Redistributions in binary form must reproduce the above copyright      #
# notice, this list of conditions and the following disclaimer in the       #
# documentation and/or other materials provided with the distribution.      #
#                                                                           #
# 3. Neither the name of the copyright holder nor the names of its          #
# contributors may be used to endorse or promote products derived from this #
# software without specific prior written permission.                       #
#                                                                           #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       #
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED #
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A           #
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER #
# OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,  #
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,       #
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR        #
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    #
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      #
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        #
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              #
#                                                                           #
#############################################################################
*/
#include "/usr/include/python3.5m/Python.h"
//#include "flagValue.h"

#include "../BioFVM/BioFVM_agent_container.h"
#include "../BioFVM/BioFVM_vector.h"
#include "../BioFVM/BioFVM.h"
#include "PhysiCell_constants.h"
#include "cell.h"
#include "sphere.h"
//I added these libraries to send *.txt files to Kafka
//#include <thread>
//#include "cppkafka/cppkafka.h"
#include "cppkafka/utils/buffered_producer.h"
#include "cppkafka/configuration.h"
//#include "../cppkafka/include/cppkafka/cppkafka.h"
//#include "cppkafka.h" //error on make
#include <iostream>
#include <fstream>

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <sys/types.h>
#include <math.h>

#define GetCurrentDir getcwd

using namespace BioFVM;
//...and this one, too.
using namespace cppkafka;
using namespace std;

//extern std::string part;

Cell_Container::Cell_Container()
{
	boundary_condition_for_pushed_out_agents = PhysiCell_constants::default_boundary_condition_for_pushed_out_agents;
	cells_ready_to_divide.clear();
	cells_ready_to_die.clear();
	num_divisions_in_current_step = 0;
	num_deaths_in_current_step = 0;
	initialized = false;
	write_all_cells = false;
	membrane_shape = "none";
	membrane_length = 0;
}	

Cell_Container::~Cell_Container()
{
}

void Cell_Container::initialize(double x_start, double x_end, double y_start, double y_end, double z_start, double z_end , double voxel_size)
{
	initialize(x_start, x_end, y_start, y_end, z_start, z_end , voxel_size, voxel_size, voxel_size);
}

void Cell_Container::initialize(double x_start, double x_end, double y_start, double y_end, double z_start, double z_end , double dx, double dy, double dz)
{
	//boundary_condition_for_pushed_out_agents = PhysiCell_constants::default_boundary_condition_for_pushed_out_agents;
	//cells_ready_to_divide.resize(0);
	//cells_ready_to_die.resize(0);

	underlying_mesh.resize(x_start, x_end, y_start, y_end, z_start, z_end , dx, dy, dz);
	agent_grid.resize(underlying_mesh.voxels.size());
	max_cell_interactive_distance_in_voxel.resize(underlying_mesh.voxels.size(), 0.0);
	agents_in_outer_voxels.resize(6);
}

int Cell_Container::nb_active_cells()
{
	int res = 0;
	for( int i=0; i < (int) all_basic_agents.size(); i++ )
	{
		if ( all_basic_agents[i]->active() )
			res ++;
	}
	return res;
}

void Cell_Container::update_all_cells(double t)
{
	update_all_cells(t, dt_settings.cell_cylce_dt_default, dt_settings.mechanics_dt_default);
}

void Cell_Container::update_all_cells(double t, double dt)
{
	update_all_cells(t, dt,dt);
}

/* Change the repulsion coefficient value of all passive cells */
void Cell_Container::set_passive_repulsion( double r )
{
		#pragma omp parallel for
		for( int i=0; i < (int) all_basic_agents.size(); i++ )
		{
			if ( all_basic_agents[i]->passive() )
				(static_cast<Sphere*>(all_basic_agents[i]))->set_repulsion(r);	
		}
}

void Cell_Container::update_all_cells_cycle( double cell_cycle_dt, double time_since_last_cycle, double t )
{
		// Reset the max_radius in each voxel. It will be filled in set_total0volume
		// It might be better if we calculate it before mechanics each time 
		std::fill(max_cell_interactive_distance_in_voxel.begin(), max_cell_interactive_distance_in_voxel.end(), 0.0);
		
		#pragma omp parallel for
		for( int i=0; i < (int) all_basic_agents.size(); i++ )
		{
			if ( !all_basic_agents[i]->passive() )
				(static_cast<Cell*>(all_basic_agents[i]))->update_cycle( cell_cycle_dt, time_since_last_cycle, t );	
		}

		// process divides / removes ( can/should be parallelized too ?) 
		#pragma omp parallel for
		for( int i=0; i < (int) cells_ready_to_divide.size(); i++ )
		{
			(cells_ready_to_divide[i])->divide();
		}
		
		for( auto cell : cells_ready_to_die )
		{	
			cell->die();	
		}

		num_divisions_in_current_step +=  cells_ready_to_divide.size();
		// The number of cells die in this time step
		num_deaths_in_current_step +=  cells_ready_to_die.size();
		
		cells_ready_to_die.clear();
		cells_ready_to_divide.clear();
}

void Cell_Container::update_all_cells_mechanics( double mechanics_dt, double time_since_last_mechanics )
{
	// Compute velocities
	#pragma omp parallel for 
	for( int i=0; i < (int) all_basic_agents.size(); i++ )
	{
		all_basic_agents[i]->update_cell_motion( time_since_last_mechanics, membrane_length, membrane_shape );
		if ( !all_basic_agents[i]->passive() )
				(static_cast<Cell*>(all_basic_agents[i]))->degrade_ecm( mechanics_dt );	
	}
		
	// Calculate new positions
	#pragma omp parallel for 
	for( int i=0; i < (int) all_basic_agents.size(); i++ )
	{
		all_basic_agents[i]->update_position(time_since_last_mechanics);
	}
		
	// Update cell indices in the container
	for( auto cell: all_basic_agents )
		cell->update_cell_index_container();
}

void Cell_Container::update_all_cells(double t, double cell_cycle_dt, double mechanics_dt)
{
	//if it is the time for running cell cycle, do it!
	double time_since_last_cycle= t- last_cell_cycle_time;

	// Means dt has to be small enough that 0.0001 is good limit ???
	if( fabs(time_since_last_cycle- cell_cycle_dt)<0.0001 || !initialized)
	{
		if(!initialized)
			time_since_last_cycle= cell_cycle_dt;
		update_all_cells_cycle( cell_cycle_dt, time_since_last_cycle, t );
		last_cell_cycle_time = t;
	}

	double time_since_last_mechanics = t - last_mechanics_time;
	
	// if( time_since_last_mechanics>= mechanics_dt || !initialized)
	if( fabs(time_since_last_mechanics - mechanics_dt)<0.0001 || !initialized)
	{
		if(!initialized)
			time_since_last_mechanics= mechanics_dt;
		update_all_cells_mechanics( mechanics_dt, time_since_last_mechanics );
		last_mechanics_time=t;
	}

	initialized=true;
	return;
}

/* Check if point is inside BM  */
int Cell_Container::inside_BM(Vector3d* pos)
{
	if ( membrane_length > 0 )
	{
		if ( membrane_shape == "duct" )
			return inside_BM_duct(pos);
		else if ( membrane_shape == "sphere" )
			return inside_BM_sphere(pos);
		else if ( membrane_shape == "sheet" )
			return inside_BM_sheet(pos);
	}
	return 0;
}

/* If point is inside BM for sphere geom */
int Cell_Container::inside_BM_sphere(Vector3d* pos)
{
	double distance_to_origin = pos->norm();  // distance to the origin 
	return ( ( membrane_length - distance_to_origin) > 0 );
}

/* If point is inside BM for sheet geom (in between two planes) */
int Cell_Container::inside_BM_sheet(Vector3d* pos)
{
	return ( ( membrane_length - fabs( (*pos)[2])) > 0 );
}

/* If point is inside BM for duct geom */
int Cell_Container::inside_BM_duct(Vector3d* pos)
{
	//Note that this function assumes that duct cap center is located at <0, 0, 0>
	if ( (*pos)[0] >= 0 ) // Cell is within the cylinder part of the duct
	{
		double distance_to_x_axis= pos->distance_to_xaxis();
		distance_to_x_axis = std::max(distance_to_x_axis, EPSILON);		// prevents division by zero
		return ( (membrane_length - distance_to_x_axis) > 0 );
	}
	
	// Cell is inside the cap of the duct
	double distance_to_origin= pos->norm();  // distance to the origin 
	return ( (membrane_length - distance_to_origin) > 0 );
}
	

/** delete the cell of given index */
void Cell_Container::delete_cell( int index )
{
	// deregister agent in from the agent container
	remove_agent( all_basic_agents[index] );
	// de-allocate (delete) the cell; 
	delete all_basic_agents[index]; 

	// performance goal: don't delete in the middle -- very expensive reallocation
	// alternative: copy last element to index position, then shrink vector by 1 at the end O(constant)
	// move last item to index location 
	if ( all_basic_agents.size() > 1 ) 
	{
		all_basic_agents[ all_basic_agents.size()-1 ]->index=index;
		all_basic_agents[index] = all_basic_agents[ all_basic_agents.size()-1 ];
		// shrink the vector
		all_basic_agents.pop_back();	
	}
	else
		all_basic_agents.resize(0); // only deleted agent was in the vector
	
}

/** add the agent to the grid */
void Cell_Container::register_agent( Basic_Agent* agent )
{
	agent_grid[agent->get_current_mechanics_voxel_index()].push_back(agent);	
}

void Cell_Container::add_new_cell( Cell* cell )
{
	all_basic_agents.push_back( cell );
	cell->index = all_basic_agents.size() - 1;
	cell->set_container(this);
}

Cell* Cell_Container::create_cell()
{
	Cell* pNew = new Cell(); 
	all_basic_agents.push_back( pNew ); 
	pNew->index = all_basic_agents.size() - 1;
	pNew->set_container(this);
	return pNew; 
}

Sphere* Cell_Container::create_sphere()
{
	Sphere* pNew = new Sphere(); 
	all_basic_agents.push_back( pNew ); 
	pNew->index = all_basic_agents.size() - 1;
	pNew->set_container(this);
	return pNew; 
}

void Cell_Container::add_agent_to_outer_voxel(Cell* agent)
{
	int escaping_face = find_escaping_face_index(agent);
	agents_in_outer_voxels[escaping_face].push_back(agent);
	agent->put_out();
}

int Cell_Container::update_voxel_in_container( double x, double y, double z )
{
	Vector3d pos(x,y,z);
	if ( is_position_valid( pos ) )
		return nearest_voxel_index( pos );
	else
		return -1;

}

bool Cell_Container::contain_any_cell(int voxel_index)
{
	return agent_grid[voxel_index].size()==0?false:true;
}

int Cell_Container::find_escaping_face_index( Cell* agent )
{
	Vector3d pos = agent->get_position();
	if (pos[0] <= underlying_mesh.bounding_box[0])
		return PhysiCell_constants::mesh_lx_face_index;
	if (pos[0] >= underlying_mesh.bounding_box[3])
		return PhysiCell_constants::mesh_ux_face_index;
	if (pos[1] <= underlying_mesh.bounding_box[1])
		return PhysiCell_constants::mesh_ly_face_index;
	if (pos[1] >= underlying_mesh.bounding_box[4])
		return PhysiCell_constants::mesh_uy_face_index;
	if (pos[2] <= underlying_mesh.bounding_box[2])
		return PhysiCell_constants::mesh_lz_face_index;
	if (pos[2] >= underlying_mesh.bounding_box[5])
		return PhysiCell_constants::mesh_uz_face_index;
	return -1;
}

void Cell_Container::flag_cell_for_division( Cell* pCell )
{ 
#pragma omp critical 
{cells_ready_to_divide.push_back( pCell );} }

void Cell_Container::flag_cell_for_removal( Cell* pCell )
{ 
#pragma omp critical 
{cells_ready_to_die.push_back( pCell );} }

int Cell_Container::writePov(double timepoint, double scale)
{
	/** \todo to move part of it in basic_agent if need pov files */
/**	std::string filename; 
	filename.resize( 1024 ); 
	sprintf( (char*) filename.c_str() , "output//cells_%i.pov" , (int)round(timepoint) ); 
	std::ofstream povFile (filename.c_str(), std::ofstream::out);
	povFile<<"#include \"colors.inc\" \n";
	povFile<<"#include \"header.inc\" \n";
	
	for( auto agent: (all_basic_agents) )
	{
		if ( !agent->passive() )
		{
		Cell* cell = static_cast<Cell*>(agent);
		std::string _nameCore;
		
		if (cell->phase_code() > 0)
		{
			int code = cell->phase_code();
			if (code ==PhysiCell_constants::Ki67_positive_premitotic || code==PhysiCell_constants::Ki67_positive_postmitotic || code==PhysiCell_constants::Ki67_positive || code==PhysiCell_constants::Ki67_negative || code==PhysiCell_constants::live)
				_nameCore="LIVE";
			else if (code==PhysiCell_constants::apoptotic)
				_nameCore="APOP";
			else if (code==PhysiCell_constants::necrotic_swelling || code==PhysiCell_constants::necrotic_lysed || code==PhysiCell_constants::necrotic)
				_nameCore="NEC";
			else if (code==PhysiCell_constants::debris)
				_nameCore="DEBR";
			else
				_nameCore="MISC";
		}
		else if( cell->type==PhysiCell_constants::TUMOR_TYPE)
			_nameCore="LIVE";
		else if( cell->type==PhysiCell_constants::VESSEL_TYPE)
			_nameCore="ENDO";
		else
			_nameCore="MISC";
	
		Vector3d pos = cell->get_position();
		std::string center= "<" + std::to_string( pos[0]/scale) + "," + std::to_string(pos[1]/scale) +","+ std::to_string(pos[2]/scale) +">";
		std::string core = "sphere {\n\t" + center + "\n\t " + std::to_string( cell->volume.radius/scale) + "\n\t FinishMacro ( " + center +","+ _nameCore+ "Finish,"+ _nameCore + "*1)\n}\n";
		povFile<< core;		
		}
	}
	
	povFile<<"#include \"footer.inc\" \n";
	povFile.close();
*/
	return 0;
}

/* The following procedure creates the cells_*.txt files */
int Cell_Container::writeCellReport(double timepoint, int cyclemode)
{
	std::string filename; 
	filename.resize( 1024 );
    std::string delimeter = ";"; // for paraview reading, doesn't like \t	
	
	// filename.c_str() <= output//cells_xxx.txt
	// !!! cells_xxx.txt has already the info!!!!!
	sprintf( (char*) filename.c_str() , "output//cells_%05d.txt" , (int)round(timepoint) ); 
	std::cout << filename.c_str() <<std::endl;

	// Output file stream povFile (=cells_xxx.txt)
	// ofstream: Stream class to write on files
	std::ofstream povFile (filename.c_str(), std::ofstream::out);
	// povFile = cells_xxx.txt

	// First line of cells_xxx.txt file
	povFile << "Time" << delimeter; 
	povFile << "ID" << delimeter;
	povFile << "x" << delimeter << "y" << delimeter << "z" << delimeter;
	povFile << "radius" << delimeter << "volume_total" << delimeter << "radius_nuclear" << delimeter << "contact_ECM" << delimeter << "freezer" << delimeter << "polarized_fraction" << delimeter << "motility" << delimeter;
    povFile	<< "cell_line" << delimeter;
	if ( cyclemode == 0 )
	   povFile << "phenotype" << delimeter << "phase" << delimeter << "elapsed_time";
	else
		povFile << "Cell_cell" << delimeter << "phase" << delimeter << "Cycle";
	//povFile << delimeter << "angle";
	povFile << delimeter << "NFkB";
	povFile << std::endl;
	//endl: end line = "\n"
	// First line of cells_xxx.txt file (the name of the attrs) have been written 

   std::ofstream passFile;
   // write passive file if option is ON
	if ( write_all_cells )
	{
		std::string passfilename; 
		passfilename.resize( 1024 );
		sprintf( (char*) passfilename.c_str() , "output//passive_cells_%05d.txt" , (int)round(timepoint) ); 
		passFile.open(passfilename.c_str(), std::ofstream::out);
		passFile << "Time" << delimeter; 
		passFile <<"ID" << delimeter;
		passFile << "x" << delimeter << "y" << delimeter << "z" << delimeter;
		passFile << "radius" << delimeter << "volume_total" << delimeter << "volume_nuclear_fluid" << delimeter << "volume_nuclear_solid" << delimeter << "volume_cytoplasmic_fluid" << delimeter << "volume_cytoplasmic_solid" << delimeter << "volume_calcified_fraction" << delimeter;
		passFile	<< "cell_line" << delimeter << "phenotype" << delimeter << "phase" << delimeter << "elapsed_time\n";
	}

	for( auto cell : all_basic_agents )
	{
		// Here we are!
		if ( !cell->passive() )
		{
			// Give the current time in cells_xxx.txt file
			povFile << timepoint << delimeter;

			// Give values to cells_xxx.txt
			cell->output( delimeter, &povFile );
		}
		else if ( write_all_cells )
		{
			passFile << timepoint << delimeter;
			cell->output( delimeter, &passFile );
		}
	}
	if ( write_all_cells )
		passFile.close();

	povFile.close(); //This is the cells_xxx.txt file with the values

	/***************** Send cells_*.txt to Kafka (the time they produced) ******************/
	// Variables declaration
	pid_t pid;
	std::string message, phase, phase_str, msg;
	int part, timestamp, normalization_value, sum_lines=-1, i=0, count_comma=0, j=0, k=0, aliveNO=0, apoptoticNO=0, necroticNO=0;
	//Open file in read mode
	std::ifstream inpovFile;
	inpovFile.open(filename.c_str());
	// Get current working directory
  	char buff[FILENAME_MAX];
  	GetCurrentDir(buff, FILENAME_MAX);
  	std::string current_working_dir(buff);
  	//cout << current_working_dir << endl;
  	int pos = current_working_dir.find("run");
  	// Take number after 'run'
  	std::string str2 = current_working_dir.substr(pos + 3);
  	// The partition where we send message
	int partno = stoi(str2);

    if(inpovFile.is_open()) {
    	cout << "\tCurrent time: " << timepoint << endl;
        //Read while data exists. Read line-by-line file-by-file (e.g. cells_00000.txt)
        while (inpovFile >> msg) {
        	// Initialize the number of commas and the phase variable
        	count_comma=0;
        	phase = ' ';
	        // Increase the number of lines
	        sum_lines++;
	        // The first line of the file (Time;ID;etc)
	        if(sum_lines == 0) {
	        	// ignore it
	        	continue;
	        }
	        // The following lines contain the info
	        //cout << msg << endl;
	        //cout << msg.length() << endl;

	        // Read the phase from the file from this line
	        for(i=0;i<msg.length();i++) {
	        	//cout << msg[i] << endl;
	        	if(msg[i]==';') {
	        		//cout << msg[i] << endl;
	        		// Count the number of delimeters we met
	        		count_comma++;
	        		//cout << "Count_comma:" << count_comma << endl;
	        		// If we reach the comma before phase (...;phase;...)
	        		if(count_comma==14) {
		        		j = i+1;
		        		k=0;
		        		//cout << "This is the first digit of the phase: " << msg[j] << endl;
		        		while(msg[j]!=';') {
		        			//cout << "Phase1: " << msg[j] << endl;
		        			phase[k] = msg[j];
		        			j++;
		        			k++;
	        			}
	        			phase_str = &phase[0];
	        			//cout << "Phase: " << phase_str << endl;
	        			// Column mapping phase
	        			if(phase_str == "0" || phase_str == "1") {
	        				aliveNO = aliveNO + 1;

	        			} else if(phase_str == "100") {
	        				apoptoticNO = apoptoticNO + 1;

	        			} else if(phase_str == "101" || phase_str == "102" || phase_str == "103") {
	        				necroticNO = necroticNO + 1;
	        			}

	        			break;
	        		}
	        	}
	        }
        }

    	// Firstly get the pid and the starting agents
    	if(timepoint == 0) {
    		// Get current PID
    		pid = getpid();
    		// Get the number of starting agents
    		normalization_value = all_basic_agents.size();
			// generate a random number in range [0,2] that defines the partition
			srand(partno);
    		part = rand() % 3;
    		cout << "Partition to send the message: " << part << endl;
    	}

        //Total amount of alive, apoptotic and necrotic cells for the current timepoint
        cout << "Alive: " << aliveNO << endl;
        cout << "Apoptotic: " << apoptoticNO << endl;
        cout << "Necrotic: " << necroticNO << endl;
        cout << "Starting Agents: " << normalization_value << endl;
        // Normalize data
        float alive_norm = (float) aliveNO/normalization_value;
        float apoptotic_norm = (float) apoptoticNO/normalization_value;
        float necrotic_norm = (float) necroticNO/normalization_value;
        // Print normalized data problem with floating point!!!!!
        cout << "Alive(normalized): " << alive_norm << endl;
        cout << "Apoptotic(normalized): " << apoptotic_norm << endl;
        cout << "Necrotic(normalized): " << necrotic_norm << endl;
        // Concatenate the pid with the current timepoint and the cells' amounts
        // Check it again.
        //if((float) timepoint==1440.02) {
        //	timestamp = (int) ceil(timepoint);
        //} else {
    	//timestamp = (int) round(timepoint);
        //}
        //Send the normalized data so as to be comparable with the data with other starting agents 
        message = std::to_string(pid) + ';' + std::to_string(timepoint) + ';' + std::to_string(alive_norm) + ';' + std::to_string(apoptotic_norm) + ';' + std::to_string(necrotic_norm);
        //cout << message << endl;
    }
    inpovFile.close();

    // message is ready to be sent to Kafka
	//Creates the topic if it doesn't exist
	MessageBuilder builder("cells1"); 
	// Define the configuration structure 
	Configuration config = { { "metadata.broker.list", "192.168.56.101:9092" } };
    // Create the producer
    BufferedProducer<string> producer(config);
    //Producer producer(config);    //Not working 
	//Produce a message
	// The message that will be sent
	cout << "Message to Kafka: " << message << endl;
	std::string str(message);
	builder.partition(0).payload(str);
    producer.add_message(builder);
    producer.flush();
    //std::this_thread::sleep_for(std::chrono::milliseconds(300));

	return 0;
}

//Creates report.txt
void Cell_Container::log_output(double t, int output_index, Microenvironment* microenvironment, std::ofstream& report_file, int cyclemode)
{
	int num_new_cells= 0;
	int num_deaths=0;
	std::cout << "current simulated time: " << t   << " minutes " << std::endl; 
	std::cout << "interval wall time: ";
	BioFVM::TOC();
	BioFVM::display_stopwatch_value( std::cout , BioFVM::stopwatch_value() ); 
	std::cout << std::endl; 
	std::cout << "total wall time: "; 
	BioFVM::RUNTIME_TOC();
	BioFVM::display_stopwatch_value( std::cout , BioFVM::runtime_stopwatch_value() ); 
	std::cout << std::endl;
	
	// Prints the current time
	std::cout << "time: "<< t <<std::endl; //t = 0, 30, 60 etc
	num_new_cells=t==0 ? all_basic_agents.size() : num_divisions_in_current_step;
	num_deaths = num_deaths_in_current_step;

	std::cout << "total number of agents (newly born, deaths): " << all_basic_agents.size() << "(" << num_new_cells << ", " << num_deaths << ")" << std::endl; 
	report_file << t << "\t" << all_basic_agents.size() << "\t" << num_new_cells << "\t" << num_deaths << "\t" << BioFVM::stopwatch_value();
   	
	int ndens = microenvironment->number_of_densities();
	for ( int i = 0; i < ndens; i++ )
	{
		report_file << "\t" << microenvironment->total_density(i);
	}
	report_file << std::endl; 
	BioFVM::TIC();
	
	num_divisions_in_current_step = 0;
	num_deaths_in_current_step = 0;
	//writePov(t, scale);
	writeCellReport(t, cyclemode);
	// Writing microenvironment state, don't for now
/**	std::string filename; 
	filename.resize( 1024 , '\0' ); 
	sprintf( (char*) filename.c_str() , "microutput/output%08d.mat" , output_index ); 
	filename.resize( strlen( filename.c_str() ) ); 
	std::cout << "\tWriting to file " << filename << " ... " << std::endl; 
	microenvironment.write_to_matlab( filename ); */
}

void Cell_Container::draw_cells_SVG( WriteSVG* svg, std::ostream* os, double zslice, std::vector<double> lims, bool plot_nuclei, int mode_color )
{
	int nagents = (int) all_basic_agents.size();
	for( int i=0; i < nagents; i++ )
	{
		if ( all_basic_agents[i]->active() )
			(static_cast<Cell*>(all_basic_agents[i]))->drawSVG( svg, os, zslice, lims, plot_nuclei, mode_color );	
	}
}
