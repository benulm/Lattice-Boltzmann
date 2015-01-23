
namespace lb
{

	void simulation::collide()
	{
//#pragma omp parallel for collapse (2)
		for (int j=0; j<static_cast<int>(l.ny); ++j)
		{
			for (int i=0; i<static_cast<int>(l.nx); ++i)
			{
				/*************** Calculating the moments ***************/
				double u = l.get_node(i,j).u();
				double v = l.get_node(i,j).v();
				double rho = l.get_node(i,j).rho();
				double u_squared = u*u + v*v;

				// M(00), M(01), M(02), M(10), M(11), ...
				std::vector<double> moments (9,0); 

				// K(0,0), K(0,1), K(0,-1), K(10), K(11), K(1,-1) ...
				std::vector<double> k_vals (9,0);

				// S(0,0), S(0,1), S(0,-1), S(10), S(11), S(1,-1) ...
				std::vector<double> s_vals (9,0); 
				// H(0,0), H(0,1), H(0,-1), H(10), H(11), H(1,-1) ...
				std::vector<double> h_vals (9,0); 

				// S(0,0), S(0,1), S(0,-1), S(10), S(11), S(1,-1) ...
				std::vector<double> s_eq_vals (9,0); 
				// H(0,0), H(0,1), H(0,-1), H(10), H(11), H(1,-1) ...
				std::vector<double> h_eq_vals (9,0); 

				// S(0,0), S(0,1), S(0,-1), S(10), S(11), S(1,-1) ...
				std::vector<double> delta_s (9,0); 
				// H(0,0), H(0,1), H(0,-1), H(10), H(11), H(1,-1) ...
				std::vector<double> delta_h_it (9,0); 

				std::vector<double> delta_h(9,0); 

				std::vector<double> f_update_it(9,0); 

				std::vector<double> f_update(9,0); 


				// for(unsigned mom=0;mom < 9; mom++){
				// 	double p;
				// 	double q;
				// 	p = (mom - (mom % 3))/3;
				// 	q = mom % 3;
				// 	double locsum = 0;
				// 	for(unsigned k=0;k<9;k++){
				// 		locsum += l.get_node(i,j).f(k)*pow(velocity_set().cx[k],p) * pow(velocity_set().cy[k],q);
				// 	}
				// 	moments[mom] = locsum / rho;
				// }

				auto n = l.get_node(i,j);
				auto c = velocity_set().c;
				auto size = velocity_set().size;
				float_type M11 = n.f(size-1)*(c[0][size-1]*c[1][size-1]);
				float_type M20 = n.f(size-1)*(c[0][size-1]*c[0][size-1]);
				float_type M02 = n.f(size-1)*(c[1][size-1]*c[1][size-1]);
				for (int i=size-2; i>=0; --i)
				{
					M11 += n.f(i)*(c[0][i]*c[1][i]);
					M20 += n.f(i)*(c[0][i]*c[0][i]);
					M02 += n.f(i)*(c[1][i]*c[1][i]);
				}
				M11 /= rho;
				M20 /= rho;
				M02 /= rho;

				// std::cout << "Bösch moments:" << std::endl;
				// std::cout << "M02:\t" << M02 << std::endl;
				// std::cout << "M20:\t" << M20 << std::endl;
				// std::cout << "M11:\t" << M11 << std::endl;

				// std::cout << "Stelli  moments:" << std::endl;
				// std::cout << "M02:\t" << moments[2] << std::endl;
				// std::cout << "M20:\t" << moments[6] << std::endl;
				// std::cout << "M11:\t" << moments[4] << std::endl;

				

				// double T  				= moments[6] + moments[2];
				// double N  				= moments[6] - moments[2];
				// double pi_xy 			= moments[4];
				// double q_xyy			= moments[5];
				// double q_yxx			= moments[7];
				// double A 				= moments[8];

				// double tilde_pi_xy 		= pi_xy - u*v;
				// double tilde_n			= N - (u*u - v*v);
				// double tilde_t			= T - (u_squared);
				// double tilde_q_xyy		= q_xyy - (2*v*tilde_pi_xy - 0.5*u*tilde_n + 0.5*u*tilde_t + u*v*v);
				// double tilde_q_yxx		= q_yxx - (2*u*tilde_pi_xy + 0.5*v*tilde_n + 0.5*v*tilde_t + u*u*v);
				// double tilde_a 			= A - (2*(u*tilde_q_xyy + v*tilde_q_yxx) + 4*u*v*tilde_pi_xy + 0.5*(u_squared)*tilde_t - 0.5*(u*u - v*v)*tilde_n + u*u*v*v);


				// double T_eq =  2*(1.0/3.0);
				// double A_eq = (1.0/9.0);

				// // Calculating s
				// for(unsigned dir=0;dir<9;dir++){
				// 	int sig = velocity_set().cx[dir];
				// 	int lam = velocity_set().cy[dir];
				// 	if(sig == 0 && lam == 0){
				// 		s_vals[dir] = rho*(4.*u*v*tilde_pi_xy - (u*u - v*v)*tilde_n/2.0 + (u_squared - 2.)*tilde_t/2.0);
				// 	}
				// 	else if(lam == 0){
				// 		s_vals[dir] = rho/2.0 * ( (1. + sig*u + u*u - v*v)/2.*tilde_n - (2*sig*v + 4.*u*v)*tilde_pi_xy + ((1.-sig*u - u_squared)/2.0)*tilde_t);
				// 	}
				// 	else if(sig == 0){
				// 		s_vals[dir] = rho/2.0 * ( (-1. - lam*v + u*u - v*v)/2.*tilde_n - (2.*lam*u + 4.*u*v)*tilde_pi_xy + ((1.-lam*v - u_squared)/2.0)*tilde_t);
				// 	}
				// 	else{
				// 		s_vals[dir] = rho / 4.0 *( (4.*u*v + sig*lam + 2.*sig*v + 2.*lam*u)*tilde_pi_xy + ((-(u*u) + v*v - sig*u + lam*v)/2.)*tilde_n + ((u_squared + sig*u + lam*v)/2.0)*tilde_t );
				// 	}
				// }
				// // Calculating s_eq
				// for(unsigned dir=0; dir<9; dir++){
				// 	int sig = velocity_set().cx[dir];
				// 	int lam = velocity_set().cy[dir];
				// 	if(sig == 0 && lam == 0){
				// 		s_eq_vals[dir] = rho*( (u_squared - 2.)*T_eq/2.0);
				// 	}
				// 	else if(lam == 0){
				// 		s_eq_vals[dir] = rho/2.0 * (  ((1-sig*u - u_squared)/2.0)*T_eq);
				// 	}
				// 	else if(sig == 0){
				// 		s_eq_vals[dir] = rho/2.0 * (  ((1-lam*v - u_squared)/2.0)*T_eq);
				// 	}
				// 	else{
				// 		s_eq_vals[dir] = rho / 4.0 *( ((u_squared + sig*u + lam*v)/2.0)*T_eq);
				// 	}
				// }

				// for(unsigned dir =0; dir<9;dir++){
				// 	delta_s[dir] = s_eq_vals[dir] - s_vals[dir];
				// }

				std::vector<double> ds (9,0);
				const float_type cPxy  = M11 - u*v;
				const float_type cN    = M20 - M02 - (u*u - v*v);
				const float_type cT    = M20 + M02 - (u*u + v*v);

				// compute delta s
				ds[0] = -(rho*(4 - (2 + 3*cN)*u*u + 24*cPxy*u*v - 2*v*v + 3*cN*v*v + 3*cT*(-2 + u*u + v*v)))/6.;
				ds[1] = (rho*(-3*cN*(1 + u + u*u - v*v) + 3*cT*(-1 + u + u*u + v*v) - 2*(-1 + u + u*u - 6*cPxy*v - 12*cPxy*u*v + v*v)))/12.;
				ds[2] = (rho*(12*cPxy*u*(1 + 2*v) + 3*cN*(1 - u*u + v + v*v) - 2*(-1 + u*u + v + v*v) + 3*cT*(-1 + u*u + v + v*v)))/12.;
				ds[3] = (rho*(-3*cN*(1 - u + u*u - v*v) + 3*cT*(-1 - u + u*u + v*v) - 2*(-1 - u + u*u + 6*cPxy*v - 12*cPxy*u*v + v*v)))/12.;
				ds[4] = (rho*(-3*cN*(-1 + u*u + v - v*v) + 3*cT*(-1 + u*u - v + v*v) - 2*(-1 + u*u + 6*cPxy*u*(1 - 2*v) - v + v*v)))/12.;
				ds[5] = -(rho*((-2 - 3*cN + 3*cT)*u + (-2 - 3*cN + 3*cT)*u*u + (-2 + 3*cN + 3*cT)*v*(1 + v) + 6*cPxy*(1 + 2*u)*(1 + 2*v)))/24.;
				ds[6] = -(rho*((2 + 3*cN - 3*cT)*u + (-2 - 3*cN + 3*cT)*u*u + (-2 + 3*cN + 3*cT)*v*(1 + v) + 6*cPxy*(-1 + 2*u)*(1 + 2*v)))/24.;
				ds[7] = -(rho*((2 + 3*cN - 3*cT)*u + (-2 - 3*cN + 3*cT)*u*u + (-2 + 3*cN + 3*cT)*(-1 + v)*v + 6*cPxy*(-1 + 2*u)*(-1 + 2*v)))/24.;
				ds[8] = -(rho*((-2 - 3*cN + 3*cT)*u + (-2 - 3*cN + 3*cT)*u*u + (-2 + 3*cN + 3*cT)*(-1 + v)*v + 6*cPxy*(1 + 2*u)*(-1 + 2*v)))/24.;


				// std::cout << "analysing difference betweeen bösch and the stellis implementation" << std::endl;
				// std::cout << "Coordinates of point: (" << i << " , " << j << ")" << std::endl;
				// std::cout <<std::setprecision(5);
				// std::cout << "Comparison of DELTA S" << std::endl;
				// std::cout << "Bösch (ds) \t Our code (dS) \t difference" << std::endl;
				// for(unsigned k=0; k<9; k++){
				// 	std::cout << ds[k] << "\t" << delta_s[k] << "\t" << ds[k]-delta_s[k] << std::endl;
				// }
				// std::cout << std::endl;

				// /******************************** Implementation Trick *******************************************/
				// Delta h_i = (f_i - f_i^{eq})- delta s_i
				std::vector<double> dh (9,0);
				dh[0] = (-2*rho*(-2 + 3*u*u + 3*v*v))/9. - n.f(0) - ds[0];
				dh[1] = (rho*(2 + 6*u + 6*u*u - 3*v*v))/18. - n.f(1) - ds[1];
				dh[2] = (rho*(2 - 3*u*u + 6*v + 6*v*v))/18. - n.f(2) - ds[2];
				dh[3] = (rho*(2 - 6*u + 6*u*u - 3*v*v))/18. - n.f(3) - ds[3];
				dh[4] = (rho*(2 - 3*u*u - 6*v + 6*v*v))/18. - n.f(4) - ds[4];
				dh[5] = (rho*(1 + 3*u*u + 3*v + 3*v*v + u*(3 + 9*v)))/36. - n.f(5) - ds[5];
				dh[6] = (rho*(1 + 3*u*u + 3*v + 3*v*v - 3*u*(1 + 3*v)))/36. - n.f(6) - ds[6];
				dh[7] = (rho*(1 + 3*u*u - 3*v + 3*v*v + u*(-3 + 9*v)))/36. - n.f(7) - ds[7];
				dh[8] = (rho*(1 + 3*u*u + u*(3 - 9*v) - 3*v + 3*v*v))/36. - n.f(8) - ds[8];

				float_type f_old[9];
				for (unsigned k =0; k<9; ++k) {
					f_old[k] = l.get_node(i,j).f(k);
				}
				velocity_set().equilibrate(l.get_node(i,j));

				// for(unsigned dir=0; dir<9; dir++){
				// 	delta_h_it[dir] = ( f_old[dir] - l.get_node(i,j).f(dir) ) - delta_s[dir];
				// }
				// std::cout <<std::setprecision(5);
				// std::cout << "Comparison of DELTA H" << std::endl;
				// std::cout << "Bösch (dh) \t Our code (dh) \t difference" << std::endl;
				// for(unsigned k=0; k<9; k++){
				// 	std::cout << dh[k] << "\t" << delta_h_it[k] << "\t" << dh[k]-delta_h_it[k] << std::endl;
				// }
				// std::cout << std::endl;
				double scpr_sh_it = 0;
				double scpr_hh_it = 0;
				for(unsigned k=0;k<9;k++){
					scpr_sh_it += (ds[k]*dh[k])/(l.get_node(i,j).f(k));
					scpr_hh_it += (dh[k]*dh[k])/(l.get_node(i,j).f(k));
				}

				double gamma_it = 1.0/beta - (2.0-1.0/beta)*(scpr_sh_it/scpr_hh_it);

				for(unsigned k=0; k<9; k++){
					l.get_node(i,j).f(k) = f_old[k] + beta*(2.0*ds[k] + gamma_it*dh[k]);
				}
				// /**************************************************************************************************/



				// /******************************* Direct implementation ********************************************/
				// // Calculating h
				// for(unsigned dir=0;dir<9;dir++){
				// 	int sig = velocity_set().cx[dir];
				// 	int lam = velocity_set().cy[dir];
				// 	if(sig == 0 && lam == 0){
				// 		h_vals[dir] = rho*(2.*u*tilde_q_xyy + 2.*v*tilde_q_yxx + tilde_a);
				// 	}
				// 	else if(lam == 0){
				// 		h_vals[dir] = rho/2.0 * (-1.0*(sig + 2.0*u)*tilde_q_xyy - 2.0*v*tilde_q_yxx - tilde_a);
				// 	}
				// 	else if(sig == 0){
				// 		h_vals[dir] = rho / 2.0 * (-1.0*(lam + 2.0*v)*tilde_q_yxx - 2.0*u*tilde_q_xyy - tilde_a);
				// 	}
				// 	else{
				// 		h_vals[dir] = rho / 4.0 *( (sig+ 2.0*u)*tilde_q_xyy + (lam + 2.0*v)*tilde_q_yxx + tilde_a);
				// 	}
				// }
				// // Calculating h_eq
				// for(unsigned dir=0;dir<9;dir++){
				// 	int sig = velocity_set().cx[dir];
				// 	int lam = velocity_set().cy[dir];
				// 	if(sig == 0 && lam == 0){
				// 		h_eq_vals[dir] = rho*(A_eq);
				// 	}
				// 	else if(lam == 0){
				// 		h_eq_vals[dir] = -1.0* rho/2.0 * (A_eq);
				// 	}
				// 	else if(sig == 0){
				// 		h_eq_vals[dir] = -1.0*rho / 2.0 * (A_eq);
				// 	}
				// 	else{
				// 		h_eq_vals[dir] = rho / 4.0 *(A_eq);
				// 	}
				// }
				// for(unsigned dir=0; dir<9; dir++){
				// 	delta_h[dir] =h_vals[dir] - h_eq_vals[dir];
				// }


				// double scpr_sh = 0;
				// double scpr_hh = 0;
				// for(unsigned k=0;k<9;k++){
				// 	scpr_sh += (delta_s[k]*delta_h[k])/(l.get_node(i,j).f(k));
				// 	scpr_hh += (delta_h[k]*delta_h[k])/(l.get_node(i,j).f(k));
				// }

				// std::vector<double> f_mirr (9,0);
				// double gamma = 1.0/beta - (2.0-1.0/beta)*(scpr_sh/scpr_hh);
				
				// for (unsigned k=0;k<9;k++){
				// 	f_mirr[k] = k_vals[k] + (2.*s_eq_vals[k] - s_vals[k]) + ((1.-gamma)*h_vals[k] + gamma*h_eq_vals[k]);
				// }

				// for (unsigned k =0; k<9; ++k){
				// 	f_update[k] = (1.0 - beta) * f_old[k] + beta*f_mirr[k];
				// }

				// /**********************************************************************************************************/

				// /************************************************* Newton Roughson for gamma ******************************/
				// // No implementation Trick
				// auto gammafunc_nit = [&](double gamma){
				// 	double ret =0;
				// 	for (unsigned k =0; k<9; ++k) {
				// 		double dh  = delta_h[k];
				// 		double ds  = delta_s[k];
				// 		double feq = l.get_node(i,j).f(k); 
				// 		ret += (dh * log(1.0+((1.0-beta*gamma)*dh - (2.0*beta-1.0)*ds)/feq));
				// 	}
				// 	return ret;
				// };
				// // Implementation Trick
				// auto gammafunc_it = [&](double gamma){
				// 	double ret =0;
				// 	for (unsigned k =0; k<9; ++k) {
				// 		double dh  = delta_h_it[k];
				// 		double ds  = delta_s[k];
				// 		double feq = l.get_node(i,j).f(k); 
				// 		ret += (dh * log(1.0+((1.0-beta*gamma)*dh - (2.0*beta-1.0)*ds)/feq));
				// 	}
				// 	return ret;
				// };

				// typedef std::pair<double, double> Result;
				// boost::uintmax_t max_iter=50;
				// boost::math::tools::eps_tolerance<double> tol(30);
				// double gamma_root_nit, gamma_root_it;
				// try{
				// 	Result r1 = boost::math::tools::toms748_solve(gammafunc_nit,0.0,5.0, tol, max_iter);
				// 	 gamma_root_nit = r1.first;
				// 	//std::cout << "alpha("<<i<<","<<j<<") = " << alpha << std::endl;
				// }
				// catch(boost::exception & e){
				// 	std::cout << "root not found -> alpha(" <<std::endl;
				// }
				// try{
				// 	Result r1 = boost::math::tools::toms748_solve(gammafunc_it,0.0,5.0, tol, max_iter);
				// 	gamma_root_it = r1.first;
				// 	//std::cout << "alpha("<<i<<","<<j<<") = " << alpha << std::endl;
				// }
				// catch(boost::exception & e){
				// 	std::cout << "root not found -> alpha(" <<std::endl;
				// }

				// /**********************************************************************************************************/

				// /**************************************** Comparison ******************************************************/

				// std::cout << "analysing difference betweeen implementation trick and \"naive\" implementation" << std::endl;
				// std::cout << "Coordinates of point: (" << i << " , " << j << ")" << std::endl;
				// std::cout <<std::setprecision(5);
				// std::cout << "Comparison of GAMMA" << std::endl;
				// std::cout << "Trick (gamma) \t Naive (gamma) \t difference" << std::endl;
				// std::cout << gamma_it << "\t\t" << gamma << "\t\t" << gamma_it-gamma << std::endl;
				// std::cout << std::endl;
				// std::cout << "Trick (gamma_root) \t Naive (gamma_root) \t difference" << std::endl;
				// std::cout << gamma_root_it << "\t\t" << gamma_root_nit << "\t\t" << gamma_root_nit-gamma_root_it << std::endl;
				// std::cout << std::endl;

				// std::cout << "Comparison of SCALAR-PRODUCT-HH" << std::endl;
				// std::cout << "Trick (gamma) \t Naive (gamma) \t difference" << std::endl;
				// std::cout << scpr_hh_it << "\t" << scpr_hh << "\t" << scpr_hh_it-scpr_hh << std::endl;
				// std::cout << std::endl;

				// std::cout << "Comparison of SCALAR-PRODUCT-SH" << std::endl;
				// std::cout << "Trick (gamma) \t Naive (gamma) \t difference" << std::endl;
				// std::cout << scpr_sh_it << "\t" << scpr_sh << "\t" << scpr_sh_it-scpr_sh << std::endl;
				// std::cout << std::endl;

				// std::cout << "Comparison of DELTA H" << std::endl;
				// std::cout << "Trick (dH) \t Naive (dH) \t difference" << std::endl;
				// for(unsigned k=0; k<9; k++){
				// 	std::cout << delta_h_it[k] << "\t" << delta_h[k] << "\t" << delta_h_it[k]-delta_h[k] << std::endl;
				// }
				// std::cout << std::endl;



				// std::cout << "Comparison of UPDATED F" << std::endl;
				// std::cout << "Trick (f') \t Naive (f') \t difference" << std::endl;
				// for(unsigned k=0; k<9; k++){
				// 	std::cout << f_update_it[k] << "\t\t" << f_update[k] << "\t\t" << f_update_it[k]-f_update[k] << std::endl;
				// }

				// std::cout << std::endl<< std::endl<< std::endl;
				
			}
		}
	}

}
