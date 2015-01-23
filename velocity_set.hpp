/** 
 *  @file
 *  @author Fabian Bösch
 *  @brief velocity set
 */

#ifndef LB_VELOCITY_SET_HPP_INCLUDED
#define LB_VELOCITY_SET_HPP_INCLUDED

#include "global.hpp"
#include <array>
#include <cmath>
#include <vector>

 namespace lb {

	struct v9;                // forward declaration
	const v9& velocity_set(); // forward declaration

	/**
	 *  @brief Lattice parameters for 9 velocity model.
	 *  
	 *  This class models a the singleton design pattern. That means there 
	 *  exists only one single instance throughout the lifetime of the 
	 *  program. To instantiate and access this object use the free function
	 *  @ref velocity_set.
	 * 
	 *  This class holds parameters like lattice weights, molecular 
	 *  velocities and speed of sound. It also exposes member functions to 
	 *  compute the equilibrium populations.
	 */
	struct v9 // singleton
	{
	private:

			/** @brief Default constructor */
		v9(){};
			/** @brief Function for instantiating the singleton is a friend */
		friend const v9& lb::velocity_set();

		float_type power(const float_type& a,const float_type& b) const{
			if (b==0.) return 1.;
			else if(b==1.) return a;
			else if(b == -1.) return 1./a;
			else std::cout << "general power operation!" << std::endl; return pow(a,b);
		}


	public:

		v9(const v9&) = delete;
		v9& operator=(const v9&) = delete;

			//                                                     0,       1,       2,       3,       4,       5,       6,       7,       8
			const std::array<float_type, 9>         W =   {{ 16.0/36,  4.0/36,  4.0/36,  4.0/36,  4.0/36,  1.0/36,  1.0/36,  1.0/36,  1.0/36}};   ///< Lattice weights

			const std::array<std::array<int, 9>, 2> c = {{{{       0,       1,       0,      -1,       0,       1,      -1,      -1,       1}}, 
				{{       0,       0,       1,       0,      -1,       1,       1,      -1,      -1}}}}; ///< Molecular velocities

				const std::array<float_type, 9> cx = {{       0,       1,       0,      -1,       0,       1,      -1,      -1,	1}};


				const std::array<float_type, 9> cy = {{       0,       0,       1,       0,      -1,       1,       1,      -1,      -1}}; 

			const float_type cs = 1.0/std::sqrt(3.0);   ///< Speed of sound

			const unsigned int size = 9;                ///< Number of velocities

			/** 
			 *  @brief Compute equilibrium.
			 * 
			 *  Compute f_eq from the locally conserved quantities rho, u and v 
			 *  (see also @ref v9::equilibrate).
			 *  @param[in,out] f_eq Pointer to an array of size 9 to store the computed values
			 *  @param[in]     rho  Local density
			 *  @param[in]     u    Local flow velocity in x-direction
			 *  @param[in]     v    Local flow velocity in y-direction
			 */

			 inline void f_eq(float_type* f_eq, float_type rho, float_type u, float_type v) const
			 {
			 	for (unsigned i = 0; i<size; ++i){
			 		f_eq[i] = rho*W[i]*(2.-sqrt(1.+3.*u*u))*(2.-sqrt(1.+3.*v*v))*power((2.*u+sqrt(1.+3.*u*u))/(1.-u),cx[i])*power((2.*v+sqrt(1.+3.*v*v))/(1.-v),cy[i]);
			 	}
				// f_eq[0] = (-2*rho*(-2 + 3*u*u + 3*v*v))/9.;
				// f_eq[1] = (rho*(2 + 6*u + 6*u*u - 3*v*v))/18.;
				// f_eq[2] = (rho*(2 - 3*u*u + 6*v + 6*v*v))/18. ;
				// f_eq[3] = (rho*(2 - 6*u + 6*u*u - 3*v*v))/18.;
				// f_eq[4] = (rho*(2 - 3*u*u - 6*v + 6*v*v))/18.;
				// f_eq[5] = (rho*(1 + 3*u*u + 3*v + 3*v*v + u*(3 + 9*v)))/36.;
				// f_eq[6] = (rho*(1 + 3*u*u + 3*v + 3*v*v - 3*u*(1 + 3*v)))/36.;
				// f_eq[7] = (rho*(1 + 3*u*u - 3*v + 3*v*v + u*(-3 + 9*v)))/36.;
				// f_eq[8] = (rho*(1 + 3*u*u + u*(3 - 9*v) - 3*v + 3*v*v))/36.;
			 }

			 inline float_type f_eq_one(float_type rho, float_type u, float_type v, int dir) const
			 {
			 	return rho*W[dir]*(2.-sqrt(1.+3.*u*u))*(2.-sqrt(1.+3.*v*v))*power((2.*u+sqrt(1.+3.*u*u))/(1.-u),cx[dir])*power((2.*v+sqrt(1.+3.*v*v))/(1.-v),cy[dir]);
			 }


			/** 
			 *  @brief Equilibrate a node.
			 * 
			 *  Compute f_eq from the locally conserved quantities rho, u and v
			 *  and set the node's population to that equilibrium ( see also 
			 *  @ref v9::f_eq).
			 *  @tparam        Node A node type
			 *  @param[in,out] n    Reference to a Node object
			 *  @param[in]     rho  Local density
			 *  @param[in]     u    Local flow velocity in x-direction
			 *  @param[in]     v    Local flow velocity in y-direction
			 */
			template <typename Node>
			 inline void equilibrate(Node& n, float_type rho, float_type u, float_type v) const
			 {
					// **************************

			 	float_type f_equil[9];
					// * fill in your code here *
			 	f_eq(f_equil, rho, u, v);
			 	for (unsigned i =0; i<size; ++i) {
			 		n.f(i) = f_equil[i];
			 	}

					// **************************
			 }

			/** 
			 *  @brief Equilibrate a node.
			 * 
			 *  Compute f_eq from the locally conserved quantities rho, u and v
			 *  and set the node's population to that equilibrium ( see also 
			 *  @ref v9::f_eq and v9::equilibrate). The locally conserved 
			 *  quantities are taken form the node object itself.
			 *  @tparam        Node A node type
			 *  @param[in,out] n    Reference to a Node object
			 */
			template <typename Node>
			 inline void equilibrate(Node& n) const
			 {
			 	return equilibrate(n, n.rho(), n.u(), n.v());
			 }
			};

	/**
	 *  @brief Get a reference single instance of the velocity set.
	 *  @return 9-velocity set
	 */
	 inline const v9& velocity_set()
	 {
	 	static v9 v_set;
	 	return v_set;
	 }

} // lb

#endif // LB_VELOCITY_SET_HPP_INCLUDED
