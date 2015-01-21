#ifndef STANDARD_COLLIDE_
#define STANDARD_COLLIDE_

namespace lb
{

	void simulation::collide()
	{
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
							ret +=
							(f_old[k]+a*(l.get_node(i,j).f(k)))*log((f_old[k]+a*(l.get_node(i,j).f(k)-f_old[k]))/velocity_set().W[k])-f_old[k]*log(f_old[k]/velocity_set().W[k]);
						}
						return ret;
					};

					typedef std::pair<double, double> Result;
					boost::uintmax_t max_iter=50000;
					boost::math::tools::eps_tolerance<double> tol(30);
					Result r1 = boost::math::tools::toms748_solve(entropy, -0.1,4., tol, max_iter);
					alpha=r1.first;
					std::cout << "alpha("<<i<<","<<j<<") = " << alpha << std::endl;
				}

				for (unsigned k =0; k<9; ++k) {
					l.get_node(i,j).f(k) = f_old[k]+alpha*beta*(l.get_node(i,j).f(k)-f_old[k]);
				}
			}
		}

	}

}

#endif // STANDARD_COLLIDE_