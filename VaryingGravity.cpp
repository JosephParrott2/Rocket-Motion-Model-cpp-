/***************
NUMERICAL MODELLING by Joseph Parrott (ID: 1413771).

Solving for the displacement and velocity of the rocket in Earth's gravitational field.
The user inputs a burntime and the code outputs the time (tau), displacement and velocity.
***************/



#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>


//Defining global constants needed for the problem. These are obtained from the worksheet or are explained in the report.
const double tg = 804.9;
const double tc = 10000;
const double g = 9.829;
const double r0 = 0.6368*pow(10,7);



//This function is the two first order ODEs obtained from the second order ODE given in the worksheet.
int func(double tau, const double y[], double f[], void *params)			
{
	double tb = *(double *)params;									//Adding the parameter for the burn time.
	f[0] = y[1];
	f[1] = tb/(tg*tg)*(tc/(1-tau) - tb/(y[0]*y[0]));
	return GSL_SUCCESS;
}



int
main (void)
{

	double A, tb, E;						//Defining variables for step size, burn time and energy respectively.
	int  s;

	FILE * fp;
	fp = fopen ("burntime3.txt", "w");

	//These allows the user to input the desired step size and the value for t_b (burn time).
	printf( "%s\n", "Please enter the desired step size" );
	std::cin >> A;
	printf( "%s\n", "Please enter enter value for burn time" );
	std::cin >> tb;

	double mu = tb;													//mu is the pointer used in the gsl algorithm and so we've set its value to our burn time																
	gsl_odeiv2_system sys = {func, 0, 2, &mu};						//inputting the function and parameter in the algorithm. The 0 is because the jacobian is not used.						

	gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk4, A, 0, 1);					//Step size A, absolute error 0, relative error 1

	double tau = 0.0;
	double y[2] = {1, 0.0};							//Setting the initial conditions (s = 1, ds/d(tau) = 0)
	


	printf("%s %36s %22s\n", "Tau", "Displacement(r_0)", "Velocity(r_0/t_b)");


//This loop actually evolves the system in time.
	while (tau <= 1)
	{	

		double E = 0.5*y[1]*y[1] - tb*tb/(tg*tg*y[0]);								//Working out the energy at each step.

		fprintf(fp, "%.15e %22.15e %22.15e\n", tau, y[0], y[1]);
		printf("%.15f %22.15f %22.15f\n", tau, y[0], y[1]);	

										
		s = gsl_odeiv2_driver_apply_fixed_step (d, &tau, A, 1, y);					//Works out the new displacement and velocity given the previous values.
		
	}

	fclose (fp);
	gsl_odeiv2_driver_free (d);

	

return 0;
}


