
#include "simulation.hpp"
#ifdef USE_OPENGL_VISUALIZATION
#include "visualization.hpp"
#endif
#include <omp.h>



int main(int argc, char *argv[])
{


	double vmax, visc, Hdivh;
	unsigned nx, h, d, w, pert;
	vmax = visc = Hdivh = nx = h = d = w = pert = 0;

	if(argc < 9){
 		std::cout << "Call the programm: " << argv[0] << "[nx, visc, v_max, h(obstacle), H/h, width(obstacle), dist(obstacle), pertubation]" << std::endl; 
 		return 0;
	}
    else {
        nx = atoi(argv[1]);
        visc = atof(argv[2]);
        vmax = atof(argv[3]);
        h = atoi(argv[4]);
        Hdivh = atoi(argv[5]);
        w = atoi(argv[6]);
        d = atoi(argv[7]);
        pert = atoi(argv[8]);
    } 

	omp_set_num_threads(std::max(omp_get_max_threads(),omp_get_num_procs()));
	
	lb::simulation* sim = new lb::simulation(nx,(unsigned) (Hdivh * h),nx*vmax/visc,vmax, visc);
	sim->initialize(w,h,d,pert);
	std::cout << *sim << std::endl;


	sim->initial_statistics_averaged();

	#ifdef USE_OPENGL_VISUALIZATION
	
		lb::visualization::initialize(sim,argc,argv);
		lb::visualization::get_instance().run();
	
	#else
	
		// Here are some hints for getting aquainted with the lattice class
		// ================================================================
		
		// how to print the lattice:
		// -------------------------
	/*	
		std::cout << sim->l << std::endl;

		// how to access the lattice:
		// --------------------------
		
		// 1) access via node proxyS
		sim->l.get_node(1,0).f(0) = 777;
		
		// 2) access data directly (make sure you know what you're doing)
		sim->l.f[0][sim->l.index(2,0)] = 3;
		
		// 3) using iterators to nodes
		(sim->l.begin() + sim->l.index(0,0))->f(0) = 1;
		
	*/
		//std::cout << sim->l << std::endl;
		
		
		// use a loop like this to run the simulation

		// std::cout << sim->l << std::endl;
	
		for (unsigned int i=0; i<1; ++i)
		{
			sim->step();
			//std::cout << "step = " << i << std::endl;

			//std::cout << sim->l << std::endl;
		}

		// std::cout << sim->l << std::endl;
		// 	for(unsigned i=0;i<sim->l.ny;i++){	
		// for(unsigned j=0;j<sim->l.nx;j++){
		
		// 		std::cout << "node (" << i << " , " << j << ") has velocities \nu: " << sim->l.get_node(i,j).u() << "\nv:" << sim->l.get_node(i,j).v() << std::endl << std::endl;
		// 	}
		// }

		//std::cout << sim->l << std::endl;
		//sim->write_vtk(2000);

	
	#endif
	
	return 0;
}
