/******************************************************************
NUMERICAL MODELLING by Joseph Parrott (ID: 1413771).

Calculating the displacement and velocities of the multistage rocket
(like previously) but loops around with many different fuel masses rather
than before when the user defined it.
******************************************************************/

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <math.h>

//Defining constants that are needed later
const double tg = 804.9;
const double tc = 10000;
const double g = 9.829;
const double U = 1000;
const double r0 = 0.6368*pow(10,7);
double m0;										//total mass
const double G = 0.667*pow(10,-10);
const double M = 0.5976*pow(10,25);

//This the first order ODEs obtained from the second order ODE in the worksheet, this is before it was put into natural units.
int func(double t, const double y[], double f[], void *params)				
{
	double r = *(double *)params;								//Here burn RATE is the parameter (r = dm/dt).
	f[0] = y[1];
	f[1] = -G*M/(y[0]*y[0]) + r*tc*g/(m0 - r*t);
	return GSL_SUCCESS;
}


int main()
{
		double r, mf2, mf1, mf, mc, mc1, mc2;					//Variables are: Burn rate, mass of fuel in first tank, mass of fuel in second tank, mass of fuel, total mass of fuel tank, mass of fuel tank one, mass of fuel tank two

		printf( "%s\n", "Please enter enter value for total mass" );
		std::cin >> m0;	
		printf( "%s\n", "Please enter enter value for rate of burn" );
		std::cin >> r;
		printf( "%s\n", "Please enter enter value for mass of fuel"); 
		std::cin >> mf;
		mc = m0 - mf;											//Total fuel tank mass = total mass - mass of fuel

		FILE * fp;
		fp = fopen ("detachvary.txt", "w");

		printf("%s %13s %32s %22s\n", "Mass of fuel in first fuel tank (kg)", "Tau", "Displacement(r_0)", "Velocity(r_0/t_b)");


		while (mf1 <= mf)						//This loops from all the fuel being in tank one, to all the fuel being in tank two
		{

		//Calculating the masses of the fuel and fuel tanks. Using fuel mass ratios to work out mc1 and mc2 (like previously)
		mf2 = mf - mf1;
		mc1 = mc * (mf1/mf);
		mc2 = mc * (mf2/mf);

		

		double mu = r;									//mu is the pointer used in the gsl algorithm and so we've set its value to our burn rate																	
		gsl_odeiv2_system sys = {func, 0, 2, &mu};		//Inputting the function and mu into the algorithm. The 0 corresponds to not using the jacobian						
						

		gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk4, 1e-5,  0, 1e-8);				//Step size 1e-5, absolute error 0, relative error 1e-8			

		double t = 0.0, T = 0.0;																						//Had to define two variables for the time, one for each loop
		double y[2] = {r0, 0.0};																						//Setting initial conditions s =1, ds/d(tau) = 0
		int i, s;

	
		//This loop is the burning of the first fuel tank.
		if (t <= mf1/r)
		{

			for (i = 0; i <= (mf1/r); i++)
			{	
		
				s = gsl_odeiv2_driver_apply_fixed_step (d, &t, 1e-3, 1000, y);	
		
			}
		}
		
			
		//This loop is the burning of the second fuel tank.	
		if (t >= mf1/r)
		{
			m0 = m0  - mf1 - mc1;										//The new mass is now the old mass - mass of the fuel that's been burnt - mass of the fuel tank that's been dropped off

			for (i = 0; i < (mf2/r); i++)
			{	
			
				s = gsl_odeiv2_driver_apply_fixed_step (d, &T, 1e-3, 1000, y);	
	
			}
				
		}
		fprintf(fp, "%.15e %22.15e %.15e %.15e\n",mf1, (t+T)/(mf/r), y[0]/r0, y[1]*(m0+mf1+mc)/(r*r0));							//Here the final displacements and velocities are plotted after the full burn is completed
		printf("%.15f %42.15f %22.15f %22.15f\n",mf1, (T+t)/(mf/r), y[0]/r0, y[1]*(m0+mf1+mc)/(r*r0));					
	
		gsl_odeiv2_driver_free (d);

		m0 = m0 + mf1 + mc1;								//Need to reset the total mass back to the original (as it was changed after the first loop)
		mf1 = mf1 + 100;									//loop around again with 100kg more fuel in the first tank

		}
	

return 0;
}
