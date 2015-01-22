
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

				// Calculating k
				for(unsigned dir=0;dir<9;dir++){
					int sig = velocity_set().cx[dir];
					int lam = velocity_set().cy[dir];
					if(sig == 0 && lam == 0){
						k_vals[dir] = rho*(1-u_squared);
					}
					else if(lam == 0){
						k_vals[dir] = rho / 2.0 * (u*u + sig*u);
					}
					else if(sig == 0){
						k_vals[dir] = rho / 2.0 * (v*v + lam*v);	
					}
					else{
						k_vals[dir] = rho / 4.0 *(sig*lam)*u*v;
					}
				}


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
				//	std::cout << "delta s(" << dir << ") is : " << delta_s[dir] << " for component (i,j) "<< i << " , " << j << std::endl; 
				}


				// Calculating h
				for(unsigned dir=0;dir<9;dir++){
					int sig = velocity_set().cx[dir];
					int lam = velocity_set().cy[dir];
					if(sig == 0 && lam == 0){
						h_vals[dir] = rho*(2*u*tilde_q_xyy + 2*v*tilde_q_yxx + tilde_a);
					}
					else if(lam == 0){
						h_vals[dir] = rho/2.0 * (-1.0*(sig + 2.0*u)*tilde_q_xyy - 2.0*v*tilde_q_yxx - tilde_a);
					}
					else if(sig == 0){
						h_vals[dir] = rho / 2.0 * (-1.0*(lam + 2.0*v)*tilde_q_yxx - 2.0*u*tilde_q_xyy - tilde_a);
					}
					else{
						h_vals[dir] = rho / 4.0 *( (sig+ 2.0*u)*tilde_q_xyy + (lam + 2.0*v)*tilde_q_yxx + tilde_a);
					}
				}
				// Calculating h_eq
				for(unsigned dir=0;dir<9;dir++){
					int sig = velocity_set().cx[dir];
					int lam = velocity_set().cy[dir];
					if(sig == 0 && lam == 0){
						h_eq_vals[dir] = rho*(A_eq);
					}
					else if(lam == 0){
						h_eq_vals[dir] = -1.0* rho/2.0 * (A_eq);
					}
					else if(sig == 0){
						h_eq_vals[dir] = -1.0*rho / 2.0 * (A_eq);
					}
					else{
						h_eq_vals[dir] = rho / 4.0 *(A_eq);
					}
				}
				for(unsigned dir=0; dir<9; dir++){
					delta_h[dir] =h_vals[dir] - h_eq_vals[dir];
				}



				float_type f_old[9];
				for (unsigned k =0; k<9; ++k) {
					f_old[k] = l.get_node(i,j).f(k);
				}
				velocity_set().equilibrate(l.get_node(i,j));

				scpr_sh = 0;
				scpr_hh = 0;
				for(unsigned k=0;k<9;k++){
					scpr_sh += (delta_s[k]*delta_h[k])/(l.get_node(i,j).f(k));
					scpr_hh += (delta_h[k]*delta_h[k])/(l.get_node(i,j).f(k));
				}

				std::vector<double> f_mirr (9,0);
				double gamma = 1.0/beta - (2.0-1.0/beta)*(scpr_sh/scpr_hh);
				
				for (unsigned k=0;k<9;k++){
					f_mirr[k] = k_vals[k] + (2.*s_eq_vals[k] - s_vals[k]) + ((1.-gamma)*h_vals[k] + gamma*h_eq_vals[k]);
				}

				for (unsigned k =0; k<9; ++k){
					l.get_node(i,j).f(k) = (1.0 - beta) * f_old[k] + beta*f_mirr[k];
				}
				
			}
		}
	}

}