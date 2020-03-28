/*****************************************
NUMERICAL MODELLING by Joseph Parrott (ID: 1413771).

Calculating the displacement and velocity of a rocket with fuel tank mass included. In this
instance the fuel is split up into two tanks and the first tank is detached after all its
fuel is burnt. The user defines the mass of the rocket, the mass of fuel and how much of 
that fuel is in the first tank. The remaining mass is the total mass of the fuel tanks and 
is split between the two tanks by taking ratios of the amount of fuel they each have. 

There are two loops: one for the first tank and one for the second. The first tank runs until
t = mf1/r because that's how long the first tank will burn for. Then the second loop goes from
t = mf1/r to t = mf2/r.

*****************************************/

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

//Defining global constants that are needed.
const double tg = 804.9;
const double tc = 10000;
const double g = 9.829;
const double U = 1000;
const double r0 = 0.6368*pow(10,7);
const double G = 0.667*pow(10,-10);
const double M = 0.5976*pow(10,25);
double m0;		


//This the first order ODEs obtained from the second order ODE in the worksheet, this is before it was put into natural units.
int func(double t, const double y[], double f[], void *params)				
{
	double r = *(double *)params;
	f[0] = y[1];
	f[1] = -G*M/(y[0]*y[0]) + r*tc*g/(m0 - r*t);
	return GSL_SUCCESS;
}

int main()
{
		double r, mf2, mf1, mf, mc,A;						//Variables are: Burn rate, mass of the second fuel tank, mass of first fuel tank, total mass of fuel tank (mf1 + mf2), step size

		printf( "%s\n", "Please enter enter value for the step size" );
		std::cin >> A;
		printf( "%s\n", "Please enter enter value for mass of rocket" );
		std::cin >> m0;	
		printf( "%s\n", "Please enter enter value for rate of burn" );
		std::cin >> r;
		printf( "%s %f\n", "Please enter enter value for mass of fuel, it must be less than ", m0 );
		std::cin >> mf;
		printf( "%s\n", "Please enter enter value for the mass of fuel in the first container (will be burnt first)" );
		std::cin >> mf1;

		mf2 = mf - mf1;										//m_(fuel in second tank) = m_(fuel) - m_(fuel in first tank)
		mc = m0 - mf;										//m_(fuel tanks) = m_(total) - m_(fuel)

		//Using the ratio of the fuel masses to work out the mass of each fuel tank
		double mc1 = mc * (mf1/mf);
		double mc2 = mc * (mf2/mf);

		printf( "%s %f\n", "The mass of container 1 is", mc1);
		printf( "%s %f\n", "The mass of container 2 is", mc2);
							
		FILE * fp;
		fp = fopen ("detach2.txt", "w");

		double mu = r;											//mu is the pointer used in the gsl algorithm and so we've set its value to our burn rate
		gsl_odeiv2_system sys = {func, 0, 2, &mu};				//Inputting the function and mu into the algorithm. The 0 corresponds to not using the jacobian						

		gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk4, A, 1e-8, 1e-8);				//Step size A, absolute error 1e-8, relative error 1e-8

		double t = 0.0, T = 0.0;								//Had to define two variables for the time, one for each loop
		double y[2] = {r0, 0.0};								//Setting initial conditions s =1, ds/d(tau) = 0
		int i, s;

		printf("%s %36s %22s\n", "Tau", "Displacement(r_0)", "Velocity(r_0/t_b)");

		
		//This loop is the burning of the first fuel tank.	
		if (t < mf1/r)
		{

			for (i = 0; i <= (mf1/r); i++)
			{	
				fprintf(fp, "%.15e %.15e %.15e\n", t/(mf/r), y[0]/r0, y[1]*mf/(r*r0));								//The units are converted to natural units here. Note that t_b = mf/r, not m0/r
				printf("%.15f %22.15f %22.15f\n", t/(mf/r), y[0]/r0, y[1]*mf/(r*r0));

				s = gsl_odeiv2_driver_apply_fixed_step (d, &t, A, 1/A, y);	
			}
		}
		
		
		//This loop is the burning of the second fuel tank.	
		if (t > mf1/r)
		{

			m0 = m0  - mf1 - mc1; 										//The new mass is now the old mass - mass of the fuel that's been burnt - mass of the fuel tank that's been dropped off.

			for (i = 0; i < (mf2/r); i++)
			{	
				fprintf(fp, "%.15e %.15e %.15e\n", (t+T)/(mf/r), y[0]/r0, y[1]*mf/(r*r0));										//Here the time is (t + T) because it has to continue from the previous loop.
				printf("%.15f %22.15f %22.15f\n", (T+t)/(mf/r), y[0]/r0, y[1]*mf/(r*r0));

				s = gsl_odeiv2_driver_apply_fixed_step (d, &T, A, 1/A, y);			
			}
		}

	fclose (fp);
	gsl_odeiv2_driver_free (d);
	

return 0;
}
