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
				 simulation(unsigned int nx, unsigned int ny, float_type Vmax_, double visc_)
				 : l(nx, ny), 
				 shift(velocity_set().size),
				 Re(nx*Vmax_/visc_), 
				 Vmax(Vmax_),
				 visc(visc_),
				 beta(1./(2.*visc/(velocity_set().cs*velocity_set().cs)+1.)),
				 time(0),
				file_output(false), // set to true if you want to write files
				output_freq(100),
				output_index(0),
				// Random number generators for pertubation
				distrx(-0.0,0.0),
				distry(-0.0,0.0),
				g1(712)
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
				}

			/**
			 *  @brief Initialize the flow field
			 *  
			 *  Initialization includes defining initial density, velocity and
			 *  populations. You can use Taylor-Green vortex flow conditions.
			 */
			 void initialize(int w, int h, int d, int perturb)
			 {

			 	coordinate<int> botwall;
			 	botwall.i=0;
			 	botwall.j=0;
			 	add_obstacle(1,l.nx,botwall);

                 add_roughness(h,w,d,perturb,perturb,perturb);

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
				/* DEBUG
			 	for (int j=0; j<static_cast<int>(l.ny); ++j)
			 	{
			 		for (int i=0; i<static_cast<int>(l.nx); ++i)
			 		{
			 			double y = l.get_node(i,j).coord.j;
			 			double x = l.get_node(i,j).coord.i;
			 			std::cout <<  "at coords(x,y) = (" << x << " , " << y << ") , the u val is: " << l.get_node(i,j).u() << std::endl;
			 			std::cout <<  "at coords(x,y) = (" << x << " , " << y << ") , the v val is: " << l.get_node(i,j).v() << std::endl;
			 		}
			 	}
				*/
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
			void wall_bc()
			{
#pragma omp parallel for
				for (unsigned int i=0; i<l.wall_nodes.size(); ++i)
				{
					for(unsigned properties=0; properties<l.wall_nodes[i].obstacle.size() ;properties++){
						int x_coord = l.wall_nodes[i].coord.i;
						int y_coord = l.wall_nodes[i].coord.j;
						switch(l.wall_nodes[i].obstacle[properties]){
							case wall_orientation::top :
							l.get_node( ((x_coord-1 + l.nx) % l.nx ), y_coord+1                   ).f(6) = l.wall_nodes[i].f(8);
							l.get_node(x_coord                      , y_coord+1                   ).f(2) = l.wall_nodes[i].f(4);
							l.get_node( ((x_coord+1) % l.nx)        , y_coord+1                   ).f(5) = l.wall_nodes[i].f(7);
							break;

							case wall_orientation::bot :
							l.get_node(((x_coord-1 + l.nx) % l.nx ), y_coord-1                   ).f(7) = l.wall_nodes[i].f(5);
							l.get_node(x_coord                     , y_coord-1                   ).f(4) = l.wall_nodes[i].f(2);
							l.get_node(((x_coord+1) % l.nx)        , y_coord-1                   ).f(8) = l.wall_nodes[i].f(6);	
							break;

							case wall_orientation::left :
							l.get_node(x_coord-1                   , ((y_coord-1 + l.ny) % l.ny)).f(7) = l.wall_nodes[i].f(5);
							l.get_node(x_coord-1                   , y_coord                    ).f(3) = l.wall_nodes[i].f(1);
							l.get_node(x_coord-1                   , ((y_coord+1)%l.ny)         ).f(6) = l.wall_nodes[i].f(8);
							break;

							case wall_orientation::right :
							l.get_node(x_coord+1                   , ((y_coord-1 + l.ny) % l.ny)).f(8) = l.wall_nodes[i].f(6);
							l.get_node(x_coord+1                   , y_coord                    ).f(1) = l.wall_nodes[i].f(3);
							l.get_node(x_coord+1                   , ((y_coord+1)%l.ny)         ).f(5) = l.wall_nodes[i].f(7);	
							break;
						}
					}
				}

			}

			void collide();

			void update_u_rho() 
			{
				float_type rho = 0;
				float_type u_temp, v_temp;

				for (int j=0; j<static_cast<int>(l.ny); ++j)
				{
					for (int i=0; i<static_cast<int>(l.nx); ++i)
					{
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

				wall_bc();
				apply_free_slip();
				
				update_u_rho();

				collide();

			 	// for (int j=0; j<static_cast<int>(l.ny); ++j)
			 	// {
			 	// 	for (int i=0; i<static_cast<int>(l.nx); ++i)
			 	// 	{
			 	// 		double y = l.get_node(i,j).coord.j;
			 	// 		double x = l.get_node(i,j).coord.i;
			 	// 			std::cout <<  "at coords(x,y) = (" << x << " , " << y << ") , the u val is: " << l.get_node(i,j).u() << std::endl;
			 	// 			std::cout <<  "at coords(x,y) = (" << x << " , " << y << ") , the v val is: " << l.get_node(i,j).v() << std::endl;
			 	// 	}
			 	// }

				

				add_body_force();

				// file io
				/*
				if ( file_output && ( ((time+1) % output_freq) == 0 || time == 0 ) )
				{
					write_fields();
					++output_index;
				}*/

				++time;

				int print = 1000;
			/*	
				if ((time-1)%print==0) {
					std::stringstream ss;
					ss << time;
					std::string ts = ss.str();
					std::string tss = "output/"+(ts)+"u_profile.dat";

					std::ofstream ofstr(tss.c_str(), std::ofstream::trunc);

					for (unsigned h=2; h<l.ny ; h++) {
						ofstr <<h <<","<<calc_stats_u_h(h) << std::endl;
					}

					ofstr.close();

				}

				int print2 = 100;
				 if ((time-1)%print2==0) {
				 	std::ofstream ofstr("output/u_h.dat", std::ofstream::app);
				 	ofstr <<time << "," << calc_stats_u_h(20) << std::endl;
				 }
				 */

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

			void calc_stats_at_h(int h);

			double calc_stats_u_h(int h);


		public: // members

			lattice l;                 ///< lattice
			std::vector<int> shift;    ///< amount of nodes to shift each population in data structure during advection
			const float_type Re;       ///< Reynolds number
			const float_type Vmax;     ///< mean flow velocity
			const float_type visc;     ///< viscosity
			const float_type beta;     ///< LB parameter beta
			unsigned int time;         ///< simulation time
			bool file_output;          ///< flag whether to write files
			unsigned int output_freq;  ///< file output frequency
			unsigned int output_index; ///< index for file naming


			double asa = -1;

			std::uniform_real_distribution<double>  distrx;
			std::uniform_real_distribution<double>  distry;
			std::mt19937 g1;
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
			float_type rho = 0;

			for (int j=0; j<static_cast<int>(l.ny); ++j)
			{
				for (int i=0; i<static_cast<int>(l.nx); ++i)
				{
					rho = l.get_node(i,j).rho();

				//double deltau = 0.568*(4.0*visc*Vmax)/(pow(l.ny,(2.0))*rho);
					double deltau = 0.5*(4.0*visc*Vmax)/(pow(l.ny,(2.0))*rho);

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

		void simulation::calc_stats_at_h(int h)
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

			std::cout << "the averaged pertubation at time t: " << time << " at height " << h << " is: I_{avg} = " << i_glob / count << std::endl;

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

} // lb


#include "standard_collide.hpp"
// #include "kbc_collide.hpp"
// #include "entropic_collide.hpp"

#endif // LB_SIMULATION_HPP_INCLUDED
