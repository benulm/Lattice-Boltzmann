
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
	
	lb::simulation* sim = new lb::simulation(nx,(unsigned) (Hdivh * h),vmax, visc,h,w,d,pert);
	//b::simulation* sim = new lb::simulation(nx,(unsigned) (Hdivh * h),vmax, visc,h);
	//lb::simulation* sim = new lb::simulation(50,50,vmax, visc);
	//sim->initialize(w,h,d,pert);
	std::cout << *sim << std::endl;


	sim->initial_statistics_averaged();

	#ifdef USE_OPENGL_VISUALIZATION
	
		lb::visualization::initialize(sim,argc,argv);
		lb::visualization::get_instance().run();
	
	#else
	


		// std::cout << sim->l << std::endl;
		std::cout << "\n\n\n\n\n" << std::endl;
		sim->step();
		// std::cout << sim->l << std::endl;
	
		// for (unsigned int i=1; i<2000001; ++i){
		// 	sim->step();
		// 	// if(i%1000 == 0)
		// 		//sim->write_vtk(i);
		// 	//std::cout << "step = " << i << std::endl;
		// }


	
	#endif
	
	return 0;
}
