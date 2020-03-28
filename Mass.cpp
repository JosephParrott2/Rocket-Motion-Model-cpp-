/***************
NUMERICAL MODELLING by Joseph Parrott (ID: 1413771).

Calculating the displacement and velocity of a rocket with fuel mass and mass of a fuel tank included.
Previously it was assumed that all the mass of the rocket is burnt. In this instance the user defines
the mass of the whole rocket and then the fraction of the total mass that is fuel. The difference between
the total mass and the fuel mass is made to be the mass of the fuel tank.

The loop runs from 0 to mf/r because now that is the burn time. But the m0 in the equation is the mass of
the fuel + the mass of the fuel tank.
***************/



#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>


//Defining global constants here as required later.
const double tg = 804.9;
const double tc = 10000;
const double g = 9.829;
const double r0 = 0.6368*pow(10,7);
const double G = 0.667*pow(10,-10);
const double M = 0.5976*pow(10,25);
double m0;									//m0 is declared but not defined as it is a user input later


//This the first order ODEs obtained from the second order ODE in the worksheet, this is before it was put into natural units.
int func(double t, const double y[], double f[], void *params)			
{
	double r = *(double *)params;							//Here burn RATE is the parameter (r = dm/dt).
	f[0] = y[1];
	f[1] = -G*M/(y[0]*y[0]) + r*tc*g/(m0 - r*t);
	return GSL_SUCCESS;
}



int
main (void)
{

	double A, r, f;									//Variables used in the user inputs below.

	FILE * fp;
	fp = fopen ("mass.txt", "w");


	//Here are the many required user inputs to work out the various masses and other needed information.
	printf( "%s\n", "Please enter the desired step size" );
	std::cin >> A;
	printf( "%s\n", "Please enter enter value for burn rate (not burn time)" );
	std::cin >> r;
	printf( "%s\n", "Please enter enter value for the total mass of the rocket" );
	std::cin >> m0;
	printf( "%s\n", "Please enter what fraction (as a decimal) of the rocket is fuel" );
	std::cin >> f;

	double mf = f*m0;											//This is the mass of the fuel in the rocket.
	double mu = r;												//Setting the value of the pointer mu as the burn rate.
															
	gsl_odeiv2_system sys = {func, 0, 2, &mu};					//Inputting variables and functions into the algorithm. Where the 0 corresponds to us not using the jacobian						

	gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk4, A, 1e-8, 1e-8);				//Step size A, absolute error 1e-8, relative error 1e-8.	

	double t = 0.0;
	double y[2] = {r0, 0.0};						//Setting initial conditions s = r0, ds/d(tau) = 0
	
	int i, s;

	printf("%s %36s %22s\n", "Tau", "Displacement(r_0)", "Velocity(r_0/t_b)");


	//This loop evolves the system from a time t = 0, to t = mf/r (all the fuel is spent).
	for (i = 0; i <= (mf/r); i++)
			{	
				//Here the units are converted back into the natural units.
				fprintf(fp, "%.15e %.15e %.15e\n", t/(mf/r), y[0]/r0, y[1]*mf/(r*r0));
				printf("%.15f %22.15f %22.15f\n", t/(mf/r),y[0]/r0, y[1]*mf/(r*r0));

				s = gsl_odeiv2_driver_apply_fixed_step (d, &t, A, 1/A, y);	
			}

	fclose (fp);
	gsl_odeiv2_driver_free (d);

	

return 0;
}


