#ifndef Anasol_
#define Anasol_
#include <fstream>
#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <iomanip>

using namespace std;


struct point
{
	point()
	: u(0), v(0), r(0)
	{}

	double u;
	double v;
	double r;
};



class Anasol
{
public:
	Anasol(int n_, double Vmax_)
		:n(n_), n2(n_*n_), lattice(n2), visc(0.001), cs(1.0/sqrt(3.0)),
		Vmax(Vmax_), beta(1./(2.*visc/(cs*cs)+1.))
	{
		const double pi(std::acos(-1.0));
		const double lambda_x 	= 1.;
		const double lambda_y 	= 1.;
		const double K_x 		= 2.*pi/n/lambda_x;
		const double K_y 		= 2.*pi/n/lambda_y;
		const double K2 		= K_x*K_x+K_y*K_y;
		const double k 			= sqrt(K2);
		const double Ma 		= Vmax/cs;

		for(unsigned j=0;j<n;j++){
			for(unsigned i=0;i<n;i++){
				operator()(i,j).u = -Vmax*K_y/k*sin(K_y*j)*cos(K_x*i);
				operator()(i,j).v = Vmax*K_x/k*sin(K_x*i)*cos(K_y*j);
				operator()(i,j).r = 1.-Ma*Ma/2./K2*(K_y*K_y*cos(2.*K_x*i)+K_x*K_x*cos(2.*K_y*j));
			}
		}
	}

	void print(int step)
	{
		const double pi(std::acos(-1.0));
		// **************************
		// * fill in your code here *
		const double lambda_x 	= 1.;
		const double lambda_y 	= 1.;
		const double K_x 		= 2.*pi/n/lambda_x;
		const double K_y 		= 2.*pi/n/lambda_y;
		const double K2 		= K_x*K_x+K_y*K_y;
		const double k 			= sqrt(K2);
		const double Ma 		= Vmax/cs;


		double t = step;

		cout << "damping after " << step << " steps is: " << exp(-visc * K2 * t) << endl;
		cout << "wirting to file" << endl;

		double damping = exp(-visc * K2 * t);
		double output_index = step;

		std::stringstream fns;
		fns << "anasol_" << std::setfill('0') << std::setw(4) << output_index << ".vtk";
		std::ofstream ofs(fns.str().c_str());
		if (ofs.is_open()){
			// Wirte header
			ofs << "# vtk DataFile Version 2.0\ntest example\nASCII\n\nDATASET STRUCTURED_POINTS\nDIMENSIONS ";
			ofs << n << " " << n << " 1" << endl;
			ofs << "ASPECT_RATIO 1 1 1\nORIGIN 0 0 0\nPOINT_DATA " << n*n << endl;
			ofs << "SCALARS point_vals double\nLOOKUP_TABLE default" << endl;

			for (unsigned int j=0; j<n; ++j){
				for (unsigned int i=0; i<n; ++i) 
				{
					ofs << damping * operator()(i,j).u << endl;
				}
			} 
		}

		cout << "file written" << endl;
	}

	point &operator() (int i, int j)
	{
		return lattice[i*n + j];
	}


private:
	int n;
	int n2;
	vector<point> lattice;


	double Vmax;
	double k2;
	double k;
	double beta;
	double visc;
	double cs;

};


#endif // Anasol_