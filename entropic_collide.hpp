#ifndef STANDARD_COLLIDE_
#define STANDARD_COLLIDE_

namespace lb
{

	void simulation::collide()
	{
#pragma omp parallel for collapse (2)
		for (int j=0; j<static_cast<int>(l.ny); ++j)
		{
			for (int i=0; i<static_cast<int>(l.nx); ++i)
			{
				float_type f_old[9];
				for (unsigned k =0; k<9; ++k) {
					f_old[k] = l.get_node(i,j).f(k);
				}


				velocity_set().equilibrate(l.get_node(i,j));

				double alpha = 2.;

				if(!l.get_node(i,j).has_flag_property("wall")){
					auto entropy = [&](double a){
						double ret =0;
						for (unsigned k =0; k<9; ++k) {
							double f_eq = l.get_node(i,j).f(k);
							double f = f_old[k];
							double ff = f+a*(f_eq-f);
							double W = velocity_set().W[k];
							ret += (ff*log(ff/W)-f*log(f/W));
						}
						return ret;
					};

					//std::cout << "(i,j)=" << i <<"\t"<< j << std::endl;
					/*
					   for (double q=0.1; q<4.1; q+=0.1){
					   std::cout << "H(" << q << ") = " << entropy(q) << std::endl;					
					   }
					 */

					double max = 0;
					for (unsigned k =0; k<9; ++k) {
						double val = fabs(f_old[k]-l.get_node(i,j).f(k))/(l.get_node(i,j).f(k));
						//std::cout << "f_old = " << f_old[k] << " f_equ = "<<l.get_node(i,j).f(k)<< " val = " << val << std::endl;
						if (val>max) max = val;
					}


					if (max>0.01) {
						typedef std::pair<double, double> Result;
						boost::uintmax_t max_iter=10;
						boost::math::tools::eps_tolerance<double> tol(30);
						try{
							Result r1 = boost::math::tools::toms748_solve(entropy, 1.,3., tol, max_iter);
							alpha = r1.first;
							//std::cout << "alpha("<<i<<","<<j<<") = " << alpha << std::endl;
						}
						catch(boost::exception & e){
							std::cout << "root not found -> alpha("<<i<<","<<j<<") = " << alpha << std::endl;
						}
					}
					else {
						//std::cout << max <<" close to equilibirium -> alpha =2" << std::endl;
						alpha =2;
					}

					for (unsigned k =0; k<9; ++k) {
						l.get_node(i,j).f(k) = f_old[k]+alpha*beta*(l.get_node(i,j).f(k)-f_old[k]);
					}
				}
			}
		}
		//std::cout << "END COLLIDE" << std::endl;		

	}
}
#endif // STANDARD_COLLIDE_
