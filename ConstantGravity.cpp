/*************************************************
NUMERICAL MODELLING by Joseph Parrott (ID: 1413771).

Solving for displacements and velocities of a rocket under
a constant gravitational field. The user inputs a burntime
and the code outputs the time (tau), displacement and velocity.
*************************************************/



#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>


//Defining global constants (actual values are discussed in report or given in the project sheet)
const double tg = 804.9;
const double tc = 10000;
const double g = 9.829;
const double U = 1;
const double r0 = 0.6368*pow(10,7);


//The system of first order ODEs derived from the second order ODE given in the work sheet
int func(double tau, const double y[], double f[], void *params)							
{
	double tb = *(double *)params;										//This is the burn time parameter which is inputted by the user
	f[0] = y[1];
	f[1] = tb*tb/(tg*tg)*(tc*U/(tb*(1-tau)) - 1);
	return GSL_SUCCESS;
}



int
main (void)
{

	double A, tb;

	FILE * fp;
	fp = fopen ("constggsl.txt", "w");

	printf( "%s\n", "Please enter the desired step size" );
	std::cin >> A;
	printf( "%s\n", "Please enter enter value for burn time" );
	std::cin >> tb;

	double mu = tb;																	
	gsl_odeiv2_system sys = {func, 0, 2, &mu};						//Here the function and burn time are inputted. The 0 is there because the jacobian is not needed								

	gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk4, A, 0, 1);					//A is the step size, 0 is the absolute error amd a relative error of 1.

	double tau = 0.0;
	double y[2] = {1, 0.0};											//Inputting initial conditions s = 1, ds/d(tau) = 0
	
	int i, s;

	printf("%s %36s %22s\n", "Tau", "Displacement(r_0)", "Velocity(r_0/t_b)");


	//This loop evolves the system from a time tau = 0 to tau = 1. 
	while (tau <= 1)
	{	
		
		fprintf(fp, "%.15e %22.15e %22.15e\n", tau, y[0], y[1]);
		printf("%.15f %22.15f %22.15f\n", tau, y[0], y[1]);	

		s = gsl_odeiv2_driver_apply_fixed_step (d, &tau, A, 1, y);					//This is what actually works out the new displacement and velocities given the time. 
		
	}

	fclose (fp);
	gsl_odeiv2_driver_free (d);

	

return 0;
}


