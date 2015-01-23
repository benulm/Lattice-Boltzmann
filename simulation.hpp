/** 
 *  @file
 *  @author Fabian Bösch
 *  @brief simulation
 */

#ifndef LB_SIMULATION_HPP_INCLUDED
#define LB_SIMULATION_HPP_INCLUDED
#include "H_root.hpp"
#include "lattice.hpp"
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <boost/math/tools/roots.hpp>

namespace lb {

	/**
	 *  @brief Simulation class implementing LB
	 * 
	 *  This class holds a lattice as member (see @ref simulation::l) and 
	 *  carries out the simulation steps on top of it. The main methods of 
	 *  this class are @ref simulation::advect() and 
	 *  @ref simulation::collide().
	 */
	class simulation
	{
		public: // ctor

			/**
			 *  @brief Construct from domain size and flow parameters
			 *  @param[in] nx    extent in x direction
			 *  @param[in] ny    extent in y direction
			 *  @param[in] _Re   Reynolds number
			 *  @param[in] _Vmax mean flow velocity
			 */
			simulation(unsigned int nx, unsigned int ny, float_type Vmax_, double visc_, int h_obstacle, int w, int d, int perturb)
				: l(nx, ny), 
				shift(velocity_set().size),
				Re(double(h_obstacle)*Vmax_/visc_), 
				Vmax(Vmax_),
				visc(visc_),
				beta(1./(2.*visc/(velocity_set().cs*velocity_set().cs)+1.)),
				time(0),
				file_output(false), // set to true if you want to write files
				output_freq(100),
				output_index(0),
				// Random number generators for pertubation
				distrx(-0.01,0.01), // 1
				distry(-0.02,0.02), // 2
				g1(712),
				// TMS boundary condition: 1 is use 0 is not use
				tms_neighbour_set(1)
		{ 
			// define amount to shift populations for advection
			for (unsigned int i=0; i<velocity_set().size; ++i)
			{
				// **************************
				// * fill in your code here *
				shift[i] = velocity_set().cx[i]*1.+velocity_set().cy[i]*l.real_nx;
				std::cout << "SHIFT:" << std::endl;
				std::cout << shift[i] << std::endl;
				// **************************
			}


			//add bottom wall
			coordinate<int> botwall;
			botwall.i=0;
			botwall.j=0;
			add_obstacle(1,l.nx,botwall);

			// Read state
			//read_state("state.dat");

			// Call initialize
			initialize();

			// Add pertubation
			add_roughness(h_obstacle,w,d,perturb,perturb,perturb);
		}


			/**
			 *  @brief Initialize the flow field
			 *  
			 *  Initialization includes defining initial density, velocity and
			 *  populations. You can use Taylor-Green vortex flow conditions.
			 */
			void initialize()
			{


				// Set initial atmospheric velocity profile : 	u = y^(1/7) + disturbance*y^(1/7)
				//												v = 0 + disturbance*y^(1/7)
				for (int j=0; j<static_cast<int>(l.ny); ++j)
				{
					for (int i=0; i<static_cast<int>(l.nx); ++i)
					{
						double y = l.get_node(i,j).coord.j;

						double vlenght = Vmax * pow((y/l.ny),1.0/7.);
						double disturb_x = distrx(g1); /* [-0.01, 0.01] */
						double disturb_y = distry(g1); /* [-0.02, 0.02] */
						l.get_node(i,j).u()   = vlenght + (vlenght * disturb_x);
						//l.get_node(i,j).u() = Vmax;
						l.get_node(i,j).v()   = 0 + (vlenght* disturb_y);
						l.get_node(i,j).rho() = 1;

						velocity_set().equilibrate(l.get_node(i,j));
					}
				}

				update_u_rho();
			} 

			/** 
			 *  @brief advect the populations
			 *  
			 *  Include periodic boundary conditions here also
			 */

			void advect()
			{
				// int forward = 0;
				// int backward = 0;


#pragma omp parallel for
				for (unsigned i=0; i<9; ++i){
					if(shift[i]>0) {
						// if(backward == 0){
						// 	std::cout << "backward loop \n\n" << std::endl;
						// }
						for (auto it = l.end()-1; it>=l.begin(); --it) {//backward loop
							// if(backward==0)
							// 	std::cout << "iterator at coords:" << it->coord.i << " , " << (*it).coord.j << std::endl;
							if(it + shift[i] < l.end())
								(it+shift[i])->f(i) = it->f(i);
						}

						// backward =1;
					}
				}


#pragma omp parallel for
				for (unsigned i=0; i<9; ++i){
					if(shift[i]<0) {
						// if(forward == 0){
						// 	std::cout << "backward loop \n\n" << std::endl;
						// }
						for (auto it = l.begin(); it<l.end(); ++it) {//forward loop
							// if(forward==0)
							// 	std::cout << "iterator at coords:" << it->coord.i << " , " << (*it).coord.j << std::endl;
							if(it + shift[i] >= l.begin())							
								(it+shift[i])->f(i) = it->f(i);
						}
						// forward = 1;
					}

				}

				//apply Boundary Conditions

#pragma omp parallel for
				for (int i = 1; i < 9; ++i) {
					if (i==1) { //east
						for (unsigned k=0; k<l.ny; ++k) {
							l.get_node(0,k).f(i) = l.get_node(l.nx,k).f(i);
						}
					}
					else if (i==2) { //north
						for (unsigned k=0; k<l.nx; ++k) {
							l.get_node(k,0).f(i) = l.get_node(k,l.ny).f(i); 
						}
					}
					else if (i==3) { //west
						for (unsigned k=0; k<l.ny; ++k) {
							l.get_node(l.nx-1,k).f(i) = l.get_node(-1,k).f(i); 
						}
					}
					else if (i==4) { //south
						for (unsigned k=0; k<l.nx; ++k) {
							l.get_node(k,l.ny-1).f(i) = l.get_node(k,-1).f(i); 
						}
					}
					else if (i==5) { //north-east
						for (unsigned k=1; k<l.nx; ++k) {
							l.get_node(k,0).f(i) = l.get_node(k,l.ny).f(i); 
						}
						for (unsigned k=1; k<l.ny; ++k) {
							l.get_node(0,k).f(i) = l.get_node(l.nx,k).f(i);
						}
						l.get_node(0,0).f(i) = l.get_node(l.nx,l.ny).f(i);
					}
					else if (i==6) { //north-west
						for (unsigned k=0; k<l.nx-1; ++k) {
							l.get_node(k,0).f(i) = l.get_node(k,l.ny).f(i); 
						}
						for (unsigned k=1; k<l.ny; ++k) {
							l.get_node(l.nx-1,k).f(i) = l.get_node(-1,k).f(i); 
						}
						l.get_node(l.nx-1,0).f(i) = l.get_node(-1,l.ny).f(i);
					}
					else if (i==7) { //south-west
						for (unsigned k=0; k<l.nx-1; ++k) {
							l.get_node(k,l.ny-1).f(i) = l.get_node(k,-1).f(i); 
						}
						for (unsigned k=0; k<l.ny-1; ++k) {
							l.get_node(l.nx-1,k).f(i) = l.get_node(-1,k).f(i); 
						}
						l.get_node(l.nx-1,l.ny-1).f(i) = l.get_node(-1,-1).f(i); 
					}
					else if (i==8) { //south-east
						for (unsigned k=1; k<l.nx; ++k) {
							l.get_node(k,l.ny-1).f(i) = l.get_node(k,-1).f(i); 
						}
						for (unsigned k=0; k<l.ny-1; ++k) {
							l.get_node(0,k).f(i) = l.get_node(l.nx,k).f(i); 
						}
						l.get_node(0,l.ny-1).f(i) = l.get_node(l.nx,-1).f(i); 					

					}

				}

			}

			/**  @brief apply wall boundary conditions */
			void bounce_back_bc()
			{
#pragma omp parallel for if(tms_neighbour_set > 1)
				for (unsigned int i=0; i<l.wall_nodes.size(); ++i)
				{
					for(unsigned properties=0; properties<l.wall_nodes[i].obstacle.size() ;properties++){
						int x_coord = l.wall_nodes[i].coord.i;
						int y_coord = l.wall_nodes[i].coord.j;

						std::vector<node*> local_wall_neighs;

						switch(l.wall_nodes[i].obstacle[properties]){
							case wall_orientation::top :
								l.get_node( ((x_coord-1 + l.nx) % l.nx ), y_coord+1                   ).f(6) = l.wall_nodes[i].f(8);								
								l.get_node(x_coord                      , y_coord+1                   ).f(2) = l.wall_nodes[i].f(4);
								l.get_node( ((x_coord+1) % l.nx)        , y_coord+1                   ).f(5) = l.wall_nodes[i].f(7);
								if(tms_neighbour_set==1){
									local_wall_neighs.push_back(&l.get_node( ((x_coord-1 + l.nx) % l.nx ), y_coord+1));
									local_wall_neighs.push_back(&l.get_node(x_coord                      , y_coord+1));
									local_wall_neighs.push_back(&l.get_node( ((x_coord+1) % l.nx)        , y_coord+1));
								}
								break;

							case wall_orientation::bot :
								l.get_node(((x_coord-1 + l.nx) % l.nx ), y_coord-1                   ).f(7) = l.wall_nodes[i].f(5);
								l.get_node(x_coord                     , y_coord-1                   ).f(4) = l.wall_nodes[i].f(2);
								l.get_node(((x_coord+1) % l.nx)        , y_coord-1                   ).f(8) = l.wall_nodes[i].f(6);	
								if(tms_neighbour_set==1){
									local_wall_neighs.push_back(&l.get_node(((x_coord-1 + l.nx) % l.nx ), y_coord-1                   ));
									local_wall_neighs.push_back(&l.get_node(x_coord                     , y_coord-1                   ));
									local_wall_neighs.push_back(&l.get_node(((x_coord+1) % l.nx)        , y_coord-1                   ));
								}
								break;

							case wall_orientation::left :
								l.get_node(x_coord-1                   , ((y_coord-1 + l.ny) % l.ny)).f(7) = l.wall_nodes[i].f(5);
								l.get_node(x_coord-1                   , y_coord                    ).f(3) = l.wall_nodes[i].f(1);
								l.get_node(x_coord-1                   , ((y_coord+1)%l.ny)         ).f(6) = l.wall_nodes[i].f(8);
								if(tms_neighbour_set==1){
									local_wall_neighs.push_back(&l.get_node(x_coord-1                   , ((y_coord-1 + l.ny) % l.ny)));
									local_wall_neighs.push_back(&l.get_node(x_coord-1                   , y_coord                    ));
									local_wall_neighs.push_back(&l.get_node(x_coord-1                   , ((y_coord+1)%l.ny)         ));
								}
								break;

							case wall_orientation::right :
								l.get_node(x_coord+1                   , ((y_coord-1 + l.ny) % l.ny)).f(8) = l.wall_nodes[i].f(6);
								l.get_node(x_coord+1                   , y_coord                    ).f(1) = l.wall_nodes[i].f(3);
								l.get_node(x_coord+1                   , ((y_coord+1)%l.ny)         ).f(5) = l.wall_nodes[i].f(7);
								if(tms_neighbour_set==1){
									local_wall_neighs.push_back(&l.get_node(x_coord+1                   , ((y_coord-1 + l.ny) % l.ny)));
									local_wall_neighs.push_back(&l.get_node(x_coord+1                   , y_coord                    ));
									local_wall_neighs.push_back(&l.get_node(x_coord+1                   , ((y_coord+1)%l.ny)         ));
								}
								break;
						}

						if(tms_neighbour_set==1){
							for (unsigned i =0; i<local_wall_neighs.size(); ++i) {
								if (!local_wall_neighs[i]->has_flag_property("wall_neigh")) {
									// set wall neighbor property
									local_wall_neighs[i]->set_flag_property("wall_neigh");
									l.wall_neigh_nodes.push_back(local_wall_neighs[i]);
								}
							}
						}

					//	std::cout << "END BOUNCE BACK" << std::endl;
						





						/*
						   node * nw_node = &l.get_node( ((x_coord-1 + l.nx) % l.nx ), ((y_coord+1)%l.ny)) ;
						   node * n_node = &l.get_node(x_coord,  ((y_coord+1)%l.ny)) ;
						   node * ne_node = &l.get_node( ((x_coord+1) % l.nx),  ((y_coord+1)%l.ny)) ;
						   node * sw_node = &l.get_node( ((x_coord-1 + l.nx) % l.nx ));
						   node * s_node = &l.get_node( x_coord, ((y_coord-1 + l.ny) % l.ny));
						   node * se_node = &l.get_node( ((x_coord+1) % l.nx),((y_coord-1 + l.ny) % l.ny));
						   node * w_node = &l.get_node( ((x_coord-1 + l.nx) % l.nx ), y_coord);
						   node * e_node = &l.get_node( ((x_coord+1) % l.nx), y_coord);

						   std::vector<node*> local_wall_neighs = {nw_node, n_node, ne_node, sw_node, s_node,se_node,w_node,e_node};
						   std::cout << "ver erstellet with size " << local_wall_neighs.size() << std::endl;	
						// check if node not yet labelled as wall neigh
						for (unsigned i =0; i<local_wall_neighs.size(); ++i) {
						if (!local_wall_neighs[i]->has_flag_property("wall_neigh"))
						{
						// set wall neighbor property
						local_wall_neighs[i]->set_flag_property("wall_neigh");
						l.wall_neigh_nodes.push_back(local_wall_neighs[i]);
						}
						//else std::cout << "already inserted! " << std::endl;
						}
						 */
						/*
						   switch(l.wall_nodes[i].obstacle[properties]){
						   case wall_orientation::top :
						   {
						   nw_node->f(6) = l.wall_nodes[i].f(8);
						   n_node->f(2) = l.wall_nodes[i].f(4);
						   ne_node->f(5) = l.wall_nodes[i].f(7);
						   break;
						   }

						   case wall_orientation::bot :
						   sw_node->f(7) = l.wall_nodes[i].f(5);
						   s_node->f(4) = l.wall_nodes[i].f(2);
						   se_node->f(8) = l.wall_nodes[i].f(6);	
						   break;

						   case wall_orientation::left :
						   sw_node->f(7) = l.wall_nodes[i].f(5);
						   w_node->f(3) = l.wall_nodes[i].f(1);
						   nw_node->f(6) = l.wall_nodes[i].f(8);
						   break;

						   case wall_orientation::right :
						   se_node->f(8) = l.wall_nodes[i].f(6);
						   e_node->f(1) = l.wall_nodes[i].f(3);
						   ne_node->f(5) = l.wall_nodes[i].f(7);	
						   break;
						   }
						 */

					}
				}

				if(tms_neighbour_set==1){
					tms_neighbour_set=2;
				}

			}

			void tms_bc() {
				//std::cout << "entering TMS " << std::endl;
				//0.) perform bounce back BC
				bounce_back_bc();

				//0.1) compute rho_target and u_target in rho and u of boundary_neigh node
				//std::cout << "Size wall neighs = " << l.wall_neigh_nodes.size() << std::endl;

				//std::cout << "Size wall = " << l.wall_nodes.size() << std::endl;
				for (unsigned i=0; i<l.wall_neigh_nodes.size(); ++i) {
					//std::cout << "wall neigh number " << i << " has coords " << l.wall_neigh_nodes[i]->coord.i << " , "<< l.wall_neigh_nodes[i]->coord.j << std::endl;
					double rho_tgt = 0;
					double u_tgt = 0;
					double v_tgt = 0;				
					
					for (unsigned k=0; k<9; ++k) {
						rho_tgt += l.wall_neigh_nodes[i]->f(k);
						u_tgt += l.wall_neigh_nodes[i]->f(k)*velocity_set().cx[k];
						v_tgt += l.wall_neigh_nodes[i]->f(k)*velocity_set().cy[k];							
					}
					l.wall_neigh_nodes[i]->rho() = rho_tgt;
					l.wall_neigh_nodes[i]->u() = u_tgt/rho_tgt;
					l.wall_neigh_nodes[i]->v() = v_tgt/rho_tgt ;
					
				}

				//std::cout << "HERE 1111111111111" << std::endl;

				//1.) set unknown velocities to f^eq(rho_tgt,u_tgt)
				for (unsigned int i=0; i<l.wall_nodes.size(); ++i)
				{
					for(unsigned properties=0; properties<l.wall_nodes[i].obstacle.size() ;properties++){
						int x_coord = l.wall_nodes[i].coord.i;
						int y_coord = l.wall_nodes[i].coord.j;

						node * nw_node = &l.get_node( ((x_coord-1 + l.nx) % l.nx ), ((y_coord+1)%l.ny)) ;
						node * n_node = &l.get_node(x_coord,  ((y_coord+1)%l.ny)) ;
						node * ne_node = &l.get_node( ((x_coord+1) % l.nx),  ((y_coord+1)%l.ny)) ;
						node * sw_node = &l.get_node( ((x_coord-1 + l.nx) % l.nx ));
						node * s_node = &l.get_node( x_coord, ((y_coord-1 + l.ny) % l.ny));
						node * se_node = &l.get_node( ((x_coord+1) % l.nx),((y_coord-1 + l.ny) % l.ny));
						node * w_node = &l.get_node( ((x_coord-1 + l.nx) % l.nx ), y_coord);
						node * e_node = &l.get_node( ((x_coord+1) % l.nx), y_coord);

						switch(l.wall_nodes[i].obstacle[properties]){
							case wall_orientation::top :
								nw_node->f(6) = velocity_set().f_eq_one(nw_node->rho(), nw_node->u(), nw_node->v(),6);
								n_node->f(2) = velocity_set().f_eq_one(n_node->rho(), n_node->u(), n_node->v(),2);
								ne_node->f(5) = velocity_set().f_eq_one(ne_node->rho(), ne_node->u(), ne_node->v(),5);
								break;

							case wall_orientation::bot :
								sw_node->f(7) = velocity_set().f_eq_one(sw_node->rho(), sw_node->u(), sw_node->v(),7);
								s_node->f(4) = velocity_set().f_eq_one(s_node->rho(), s_node->u(), s_node->v(),4);
								se_node->f(8) = velocity_set().f_eq_one(se_node->rho(), se_node->u(), se_node->v(),8);
								break;

							case wall_orientation::left :
								sw_node->f(7) = velocity_set().f_eq_one(sw_node->rho(), sw_node->u(), sw_node->v(),7);
								w_node->f(3) = velocity_set().f_eq_one(w_node->rho(), w_node->u(), w_node->v(),3);
								nw_node->f(6) = velocity_set().f_eq_one(nw_node->rho(), nw_node->u(), nw_node->v(),6);
								break;

							case wall_orientation::right :
								se_node->f(8) = velocity_set().f_eq_one(se_node->rho(), se_node->u(), se_node->v(),8);
								e_node->f(1) = velocity_set().f_eq_one(e_node->rho(), e_node->u(), e_node->v(),1);
								ne_node->f(5) = velocity_set().f_eq_one(ne_node->rho(), ne_node->u(), ne_node->v(),5);
								break;
						}

					}
				}

				//1.1) forall i: f_i = f_i + f_i^eq(rho_tgt, u_tgt)-f_eq(rho_loc,u_loc)
				for (unsigned i=0; i<l.wall_neigh_nodes.size(); ++i) {
					//store rho_loc and u_loc 
					double rho_loc = 0;
					double u_loc = 0;
					double v_loc = 0;
					for (unsigned k=0; k<9; ++k) {
						rho_loc += l.wall_neigh_nodes[i]->f(k);
						u_loc += l.wall_neigh_nodes[i]->f(k)*velocity_set().cx[k];
						v_loc += l.wall_neigh_nodes[i]->f(k)*velocity_set().cy[k];
					}
					u_loc /= rho_loc;
					v_loc /= rho_loc;
					//2.)
					float_type f_eq_tgt[9];
					float_type f_eq_loc[9];

					// compute f_eq_local and f_eq_target
					velocity_set().f_eq(f_eq_tgt,l.wall_neigh_nodes[i]->rho(),	l.wall_neigh_nodes[i]->u(), l.wall_neigh_nodes[i]->v());
					velocity_set().f_eq(f_eq_loc,rho_loc,u_loc,v_loc);

					for (unsigned k=0; k<9; ++k) {
						l.wall_neigh_nodes[i]->f(k)+= (f_eq_tgt[k]-f_eq_loc[k]);
					}
				}

			}

			void collide();

			void update_u_rho() 
			{

#pragma omp parallel for collapse (2)
				for (int j=0; j<static_cast<int>(l.ny); ++j)
				{
					for (int i=0; i<static_cast<int>(l.nx); ++i)
					{
						float_type rho = 0;
						float_type u_temp, v_temp;

						//update rho
						rho = 0;
						u_temp = 0; v_temp = 0;

						for (unsigned k=0; k<9; ++k) {
							rho += l.get_node(i,j).f(k);
						}
						l.get_node(i,j).rho() = rho;

						//update velocity
						for (unsigned k=0; k<velocity_set().size; ++k) {
							u_temp += l.get_node(i,j).f(k)*velocity_set().cx[k];
							v_temp += l.get_node(i,j).f(k)*velocity_set().cy[k];
						}

						u_temp /= rho;
						v_temp /= rho;

						l.get_node(i,j).u() = u_temp;
						l.get_node(i,j).v() = v_temp;

					}
				}

			}

			/** @brief LB step */
			void step()
			{
				//std::cout << "STEP starts" << std::endl;	
				advect();
				if(tms_neighbour_set == 0){
					bounce_back_bc();
				}
				else{
					tms_bc();
				}
				apply_free_slip();

				update_u_rho();

				collide();

				add_body_force();

				++time;
				int print = 50000;

				if ((time-1)%print==0){
					print_u_profile(time);
					print_I_profile(time);
					print_u_h(time, l.nx/2);
					print_state(time);
				}

			}


			void add_obstacle(int height, int width, coordinate<int> anker)
			{
				int x_start = anker.i;
				int y_start = anker.j;
				bool top = false;
				if(y_start==0){
					top = true;
				}

				for(int j=0; j<height; j++){
					for(int i=0; i<width; i++){
						l.get_node(x_start+i, y_start+j).obstacle.clear();
						if(top){
							l.get_node(x_start+i, y_start+j).obstacle.push_back(wall_orientation::top);
						}
						else{
							l.get_node(x_start+i, y_start+j).obstacle.push_back(wall_orientation::bot);
						}
						if(i==0){
							l.get_node(x_start+i, y_start+j).obstacle.push_back(wall_orientation::left);
						}
						else if(i == (width-1)){
							l.get_node(x_start+i, y_start+j).obstacle.push_back(wall_orientation::right);
						}
					}
				}

				coordinate<int> maxcoord;
				maxcoord.i = x_start + width-1;
				maxcoord.j = y_start + height-1;

				l.add_wall(anker, maxcoord);
			}

		public: // write to file

			/** write macroscopic variables to ascii file */
			void write_fields()
			{
				std::stringstream fns;
				fns << "output/data_" << std::setfill('0') << std::setw(4) << output_index << ".txt";
				l.write_fields(fns.str());
			}

			/** write vtk file */
			void write_vtk(int step)
			{
				//update_u_rho();
				std::stringstream fns;
				fns << "output/data_" << std::setfill('0') << std::setw(4) << step << ".vtk";
				l.write_vtk(fns.str());
			}


		public: // print

			/** print to output stream */
			friend std::ostream& operator<<(std::ostream& os, const simulation& sim)
			{
				os << "simulation parameters\n" 
					<< "---------------------\n";
				os << "domain: " << sim.l.nx << " x " << sim.l.ny << "\n";
				os << "Re:     " << sim.Re << "\n";
				os << "Vmax:   " << sim.Vmax << "\n";
				os << "visc:   " << sim.visc << "\n";
				os << "beta:   " << sim.beta << "\n";
				return os;
			}


		public: // Boundary Conditions

			void apply_free_slip();


		public: // Force calculations

			void add_body_force();


		public: // Obstacles

			void add_roughness(int height, int width, int dist, int width_tol, int height_tol, int dist_tol);


		public: // Statistics

			void initial_statistics_averaged();

			void calc_stats_averaged();

			int average_speed_averaged();

			double calc_stats_at_h(int h);

			double calc_stats_u_h(int h);

		public: //Output
			void print_u_profile(int t);

			void print_I_profile(int t);

			void print_u_h(int t, int h);

			void print_state(int t);

			void read_state(std::string s);



		public: // members

			lattice l;                 ///< lattice
			std::vector<int> shift;    ///< amount of nodes to shift each population in data structure during advection
			float_type Re;       		///< Reynolds number
			float_type Vmax;     		///< mean flow velocity
			float_type visc;    		///< viscosity
			float_type beta;     		///< LB parameter beta
			unsigned int time;         ///< simulation time
			bool file_output;          ///< flag whether to write files
			unsigned int output_freq;  ///< file output frequency
			unsigned int output_index; ///< index for file naming


			double asa = -1;

			std::uniform_real_distribution<double>  distrx;
			std::uniform_real_distribution<double>  distry;
			std::mt19937 g1;

			// Flag for using TMS boundary-conditions
			int tms_neighbour_set;	
	};





	/////////////////////////////////////////////////////////////////////////
	// Force term:
	// delta u = F/rho dt							(1)
	// (half) pipe flow: dp/dz= -4*u_max*mu/(h^2)	(2)
	// pressure gradient force: a = -1/rho*dp/dz	(3)
	// (1),(2) in (3): delta u = 4*nu/rho*u_max/(h^2)*dt
	/////////////////////////////////////////////////////////////////////////
	void simulation::add_body_force()
	{
#pragma omp parallel for collapse (2)
		for (int j=0; j<static_cast<int>(l.ny); ++j)
		{
			for (int i=0; i<static_cast<int>(l.nx); ++i)
			{

				float_type rho = 0;
				rho = l.get_node(i,j).rho();

				//double deltau = 0.568*(4.0*visc*Vmax)/(pow(l.ny,(2.0))*rho);
				double deltau = (4.0*visc*Vmax)/(pow(l.ny,(2.0)));

				//double deltau = (4.0*visc*Vmax)/(pow(l.ny,(2.0))*rho);
				//std::cout << "Delta u = " << deltau << std::endl;
				//double deltau = 1e-6;

				float_type f_old[9];
				for (unsigned k =0; k<9; ++k) {
					f_old[k] = l.get_node(i,j).f(k);
				}

				velocity_set().equilibrate(l.get_node(i,j),rho,l.get_node(i,j).u() + deltau, l.get_node(i,j).v());

				for (unsigned k =0; k<9; ++k) {
					f_old[k] += l.get_node(i,j).f(k);
				}

				velocity_set().equilibrate(l.get_node(i,j),rho,l.get_node(i,j).u(), l.get_node(i,j).v());

				for (unsigned k =0; k<9; ++k) {
					l.get_node(i,j).f(k) = f_old[k] - l.get_node(i,j).f(k);
				}				
			}
		}
	}

	//surface roughness
	void simulation::add_roughness(int height, int width, int dist, int width_tol = 0, int height_tol = 0, int dist_tol	= 0)
	{
		std::default_random_engine generator(712);

		std::uniform_int_distribution<int> width_perturb(-width_tol,width_tol);
		std::uniform_int_distribution<int> height_perturb(-height_tol,height_tol);
		std::uniform_int_distribution<int> dist_perturb(-dist_tol,dist_tol);			

		int pwidth, pheight, pdist;

		//calculate the perturbed quantities
		pwidth = width+width_perturb(generator);
		pheight = height+height_perturb(generator);
		pdist = dist+dist_perturb(generator);

		coordinate<int> anchor;
		anchor.i = pdist;
		anchor.j = 0;

		while (anchor.i + pwidth <= (int)l.nx) {
			add_obstacle(pheight, pwidth, anchor);
			anchor.i = anchor.i + pwidth + pdist;

			pwidth = width+width_perturb(generator);
			pheight = height+height_perturb(generator);
			pdist = dist+dist_perturb(generator);
		}
	}


	//calculating statistics:
	int simulation::average_speed_averaged()
	{
		if(asa==-1){
			double total_speed = 0;
			int count = 0;
			for(unsigned i=0; i<l.ny; i++){
				total_speed += Vmax * pow((i/l.ny),1.0/7);
				count++;
			}
			total_speed /= count;

			asa = total_speed;

		}
		return asa;

	}

	void simulation::initial_statistics_averaged()
	{
		if(time !=0){
			std::cout << "time already stepped, no initial average possible" << std::endl;
			return;
		}
		else{
			double i_glob = 0;
			int count = 0;
			for(unsigned i=0; i<l.ny; i++){
				if(!l.get_node(0,i).has_flag_property("wall")){
					double u_initial = Vmax * pow( ((double)i)/l.ny,1.0/7);
					double i_loc = (sqrt(pow(l.get_node(0,i).v(),2) + pow(l.get_node(0,i).u() - u_initial,2)) )/ u_initial;
					i_glob += i_loc;
					count++;
				}
			}

			std::cout << "the averaged initial pertubation is: I_{avg} = " << i_glob / count << std::endl;
		}
	}

	void simulation::calc_stats_averaged()
	{
		double i_glob = 0;
		int count = 0;
		for(unsigned i=0; i<l.ny; i++){
			if(!l.get_node(0,i).has_flag_property("wall")){
				double u_initial = Vmax * pow( ((double)i)/l.ny,1.0/7);
				double i_loc = (sqrt(pow(l.get_node(0,i).v(),2) + pow(l.get_node(0,i).u() - u_initial,2)) )/ u_initial;
				i_glob += i_loc;
				count++;
			}
		}

		std::cout << "the averaged pertubation at time t: " << time << " is: I_{avg} = " << i_glob / count << std::endl;

	}

	double simulation::calc_stats_at_h(int h)
	{
		double i_glob = 0;
		int count = 0;
		for(unsigned i=0; i<l.nx; i++){
			if(!l.get_node(i,h).has_flag_property("wall")){
				double u_initial = Vmax * pow( ((double)h)/l.ny,1.0/7);
				double i_loc = (sqrt(pow(l.get_node(i,h).v(),2) + pow(l.get_node(i,h).u() - u_initial,2)) )/ u_initial;
				//std::cout << "local i is: " << i_loc << std::endl;
				i_glob += i_loc;
				count++;
			}
		}
		return i_glob/count;
		//std::cout << "the averaged pertubation at time t: " << time << " at height " << h << " is: I_{avg} = " << i_glob / count << std::endl;

	}

	double simulation::calc_stats_u_h(int h)
	{
		double u_glob = 0;
		int count = 0;
		for(unsigned i=0; i<l.nx; i++){
			if(!l.get_node(i,h).has_flag_property("wall")){
				u_glob += l.get_node(i,h).u();
				count++;
			}
		}

		return u_glob/count;

	}


	void simulation::apply_free_slip()
	{
		for(unsigned i=0;i<l.nx;i++){
			l.get_node(i,l.ny-1).f(4) = l.get_node(i,l.ny).f(2);
			l.get_node(i,l.ny-1).f(8) = l.get_node(i+1,l.ny).f(5);
			l.get_node(i,l.ny-1).f(7) = l.get_node(i-1,l.ny).f(6);
		}
	}


	//Output functions:
	void simulation::print_u_profile(int t){
		std::stringstream ss;
		ss << t;
		std::string ts = ss.str();
		std::string tss = "output/"+(ts)+"u_profile.dat";

		std::ofstream ofstr(tss.c_str(), std::ofstream::trunc);

		for (unsigned h=2; h<l.ny ; h++) {
			ofstr <<h <<","<<calc_stats_u_h(h) << std::endl;
		}

		ofstr.close();

	}


	void simulation::print_I_profile(int t){				
		std::stringstream ss;
		ss << t;
		std::string ts = ss.str();
		std::string tss = "output/"+(ts)+"I_profile.dat";
		std::ofstream ofstr(tss.c_str(), std::ofstream::trunc);

		for (unsigned h=2; h<l.ny ; h++) {
			ofstr <<h <<","<<calc_stats_at_h(h) << std::endl;
		}
		ofstr.close();

	}



	void simulation::print_u_h(int t, int h) 
	{
		std::ofstream ofstr("output/u_h.dat", std::ofstream::app);
		ofstr <<t << "," << calc_stats_u_h(h) << std::endl;
	}

	void simulation::print_state(int t) {
		std::stringstream ss;
		ss << t;
		std::string ts = ss.str();
		std::string tss = "output/"+(ts)+"state.dat";
		std::ofstream ofstr(tss.c_str(), std::ofstream::trunc);
		ofstr << "State of simulation to be used as initial state: \n";
		ofstr << "V_max = " << Vmax << std::endl;
		ofstr << "nu = " << visc << std::endl;
		ofstr << "system size = [ " << l.nx << " , " << l.ny << " ]" << std::endl;
		ofstr << "################################################################\n";
		ofstr << "i\tj\tf0\tf1\tf2\tf3\tf4\tf5\tf6\tf7\tf8\tu\tv\trho\n";

		for (int j=0; j<static_cast<int>(l.ny); ++j)
		{
			for (int i=0; i<static_cast<int>(l.nx); ++i)
			{
				ofstr << i << "\t" << j << "\t";
				for (unsigned k=0; k<9; ++k) {
					ofstr << l.get_node(i,j).f(k) << "\t";
				}
				ofstr << l.get_node(i,j).u() << "\t" << l.get_node(i,j).v() << "\t" << l.get_node(i,j).rho() << std::endl;
			}
		}

		ofstr.close();
	}

	void simulation::read_state(std::string s)
	{
		std::ifstream simfile (s.c_str());
		char title[300];
		std::string buffer;
		double vmax, visc,u,v,rho;
		std::vector<double> fvals (9,0);
		unsigned lnx,lny, i_val, j_val;
		simfile.getline (title,300);
		simfile >> buffer >> buffer >> vmax;
		std::cout << buffer << std::endl;
		simfile >> buffer >> buffer >> visc;
		std::cout << buffer << std::endl;
		simfile >> buffer >> buffer >> buffer >> buffer >> lnx >> buffer >> lny >> buffer;
		std::cout << buffer << std::endl;
		simfile >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer;
		this->Vmax = vmax;
		this->visc = visc;
		this->beta = 1./(2.*visc/(velocity_set().cs*velocity_set().cs)+1.);
		if(lnx != l.nx || lny != l.ny){
			std::cout << "error: size missmatch, deisred size was (" << lnx <<" , " << lny << std::endl;
			return;
		}
		std::cout << buffer << std::endl;

		for (unsigned j=0; j<lny; ++j)
		{
			for (unsigned i=0; i<lnx; ++i)
			{
				simfile >> i_val >> j_val;

				for (unsigned k=0; k<9; ++k) {
					simfile >> fvals[k];
				}
				simfile >> u;
				simfile >> v;
				simfile >> rho;
				for(unsigned k=0;k<9;k++){
					l.get_node(i_val,j_val).f(k) = fvals[k]; 
				}
				l.get_node(i_val,j_val).u() = u;
				l.get_node(i_val,j_val).v() = v;
				l.get_node(i_val,j_val).rho() = rho;
			}
		}

	}



} // lb


//#include "standard_collide.hpp"
//#include "kbc_collide.hpp"
#include "entropic_collide.hpp"

#endif // LB_SIMULATION_HPP_INCLUDED
