
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

				// S(0,0), S(0,1), S(0,-1), S(10), S(11), S(1,-1) ...
				std::vector<double> s_vals (9,0); 
				// S(0,0), S(0,1), S(0,-1), S(10), S(11), S(1,-1) ...
				std::vector<double> s_eq_vals (9,0); 

				// S(0,0), S(0,1), S(0,-1), S(10), S(11), S(1,-1) ...
				std::vector<double> delta_s (9,0); 
				// H(0,0), H(0,1), H(0,-1), H(10), H(11), H(1,-1) ...
				std::vector<double> delta_h (9,0); 


				for(unsigned mom=0;mom < 9; mom++){
					double p;
					double q;
					p = (mom - (mom % 3))/3;
					q = mom % 3;
					double locsum = 0;
					for(unsigned k=0;k<9;k++){
						locsum += l.get_node(i,j).f(k)*pow(velocity_set().cx[k],p) * pow(velocity_set().cy[k],q);
					}
					moments[mom] = locsum / rho;
				}
				

				double T  				= moments[6] + moments[2];
				double N  				= moments[6] - moments[2];
				double pi_xy 			= moments[4];
				double q_xyy			= moments[5];
				double q_yxx			= moments[7];
				double A 				= moments[8];

				double tilde_pi_xy 		= pi_xy - u*v;
				double tilde_n			= N - (u*u - v*v);
				double tilde_t			= T - (u_squared);
				double tilde_q_xyy		= q_xyy - (2*v*tilde_pi_xy - 0.5*u*tilde_n + 0.5*u*tilde_t + u*v*v);
				double tilde_q_yxx		= q_yxx - (2*u*tilde_pi_xy + 0.5*v*tilde_n + 0.5*v*tilde_t + u*u*v);
				double tilde_a 			= A - (2*(u*tilde_q_xyy + v*tilde_q_yxx) + 4*u*v*tilde_pi_xy + 0.5*(u_squared)*tilde_t - 0.5*(u*u - v*v)*tilde_n + u*u*v*v);


				double T_eq =  2*(1.0/3.0);
				double A_eq = (1.0/9.0);

				// Calculating s
				for(unsigned dir=0;dir<9;dir++){
					int sig = velocity_set().cx[dir];
					int lam = velocity_set().cy[dir];
					if(sig == 0 && lam == 0){
						s_vals[dir] = rho*(4*u*v*tilde_pi_xy - (u*u - v*v)*tilde_n/2.0 + (u_squared - 2)*tilde_t/2.0);
					}
					else if(lam == 0){
						s_vals[dir] = rho/2.0 * ( (1 + sig*u + u*u - v*v)/2*tilde_n - (2*sig*v + 4*u*v)*tilde_pi_xy + ((1-sig*u - u_squared)/2.0)*tilde_t);
					}
					else if(sig == 0){
						s_vals[dir] = rho/2.0 * ( (-1 - lam*v + u*u - v*v)/2*tilde_n - (2*lam*u + 4*u*v)*tilde_pi_xy + ((1-lam*v - u_squared)/2.0)*tilde_t);
					}
					else{
						s_vals[dir] = rho / 4.0 *( (4*u*v + sig*lam + 2*sig*v + 2*lam*u)*tilde_pi_xy + ((-(u*u) + v*v - sig*u + lam*v)/2)*tilde_n + ((u_squared + sig*u + lam*v)/2.0)*tilde_t );
					}
				}


				double T_eq =  2*(1.0/3.0);
				// Calculating s_eq
				for(unsigned dir=0; dir<9; dir++){
					int sig = velocity_set().cx[dir];
					int lam = velocity_set().cy[dir];
					if(sig == 0 && lam == 0){
						s_eq_vals[dir] = rho*( (u_squared - 2)*T_eq/2.0);
					}
					else if(lam == 0){
						s_eq_vals[dir] = rho/2.0 * (  ((1-sig*u - u_squared)/2.0)*T_eq);
					}
					else if(sig == 0){
						s_eq_vals[dir] = rho/2.0 * (  ((1-lam*v - u_squared)/2.0)*T_eq);
					}
					else{
						s_eq_vals[dir] = rho / 4.0 *( ((u_squared + sig*u + lam*v)/2.0)*T_eq);
					}
				}

				for(unsigned dir =0; dir<9;dir++){
					delta_s[dir] = s_vals[dir] - s_eq_vals[dir] ;
				}

				// Delta h_i = (f_i - f_i^{eq})- delta s_i
				float_type f_old[9];
				for (unsigned k =0; k<9; ++k) {
					f_old[k] = l.get_node(i,j).f(k);
				}
				velocity_set().equilibrate(l.get_node(i,j));

				for(unsigned dir=0; dir<9; dir++){
					delta_h[dir] = ( f_old[dir] - l.get_node(i,j).f(dir) ) - delta_s[dir];
				}

				double scpr_sh = 0;
				double scpr_hh = 0;
				for(unsigned k=0;k<9;k++){
					scpr_sh += (delta_s[k]*delta_h[k])/(l.get_node(i,j).f(k));
					scpr_hh += (delta_h[k]*delta_h[k])/(l.get_node(i,j).f(k));
				}

				double gamma = 1.0/beta - (2.0-1.0/beta)*(scpr_sh/scpr_hh);

				for(unsigned k=0; k<9; k++){
					l.get_node(i,j).f(k) = f_old[k] - beta*(2.0*delta_s[k] + gamma*delta_h[k]);
				}
				
			}
		}
	}

}