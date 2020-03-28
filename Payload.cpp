/***************
NUMERICAL MODELLING by Joseph Parrott (ID: 1413771).

This code calculates the energy of the system at each intervel, if it remains negative the process goes on and if it becomes positive, the values are printed out and the loop is terminated. This is repeated for t_b values between 100 and t_c in steps of 100.

Note that the energy is plotted out. Ideally these need to be as close to 0 as possible but showing the values
shows the user the slight inaccuracy due it taking too large of a step.
***************/



#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>


//Some variables are declared globally here
double tg;
double tc;
double g ;
const double r0 = 0.6368*pow(10,7);


//This function contains the system of first order ODEs obtained from the second order ODE in the worksheet
int func(double tau, const double y[], double f[], void *params)				
{
	double tb = *(double *)params;							//Using burn time as a parameter
	f[0] = y[1];
	f[1] = tb/(tg*tg)*(tc/(1-tau) - tb/(y[0]*y[0]));
	return GSL_SUCCESS;
}



int
main (void)
{

	double P;
	//Allowing the user to choose if they wish to work in the context of Jupiter or Earth
	printf("%s\n", "Would you like to work out the payload for Jupiter or Earth?");
	printf("%s\n", "[1] Earth");
	printf("%s\n", "[2] Jupiter");
	std::cin >> P;
	FILE * fp;

	//Defining the variables corresponding to Earth
	if (P == 1)
	{
		fp = fopen ("earthpayload.txt", "w");
		printf("%s\n", "You have selected Earth");
		tg = 804.9;
		tc = 10000;
		g = 9.829;
	}

	//Defining the variables corresponding to Jupiter
	if (P == 2)
	{
		fp = fopen ("jupiterpayload.txt", "w");
		printf("%s\n", "You have selected Jupiter");
		tc = 100000/24.79;
		g = 24.79;
		tg = 1679;
	}

	double A, tb = 100;
	printf( "%s\n", "Please enter the desired step size, recommend something smaller than 1e-4 for good results" );
	std::cin >> A;

	printf("%s %18s %35s %9s\n", "t_b(s)", "Tau", "Energy((r_0/t_b)^2)", "Payload");


	//This loop repeats the whole process from different t_b values.
	while (tb <= tc)
	{

		int s;
		double mu = tb;						//Setting the pointer mu to the value of t_b
																	
		gsl_odeiv2_system sys = {func, 0, 2, &mu};					//Inputting values into the algorithm. With the 0 corresponding to there being no jacobian.		

		gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk4, A, 0, 1);				//Step size A, absolute error 0, relative error 1.

		double tau = 0.0;
		double y[2] = {1, 0.0};							//Setting the initial conditions s = 1, ds/d(tau) = 0.
	
		

	
		while (tau <= 1)
		{	
			double E = 0.5*y[1]*y[1] - tb*tb/(tg*tg*y[0]);				//Works out the energy as required by the if statements

			if (E < 0)													//When energy is negative (i.e. escape velocity is not achieved) continue the process
			{
				
				s = gsl_odeiv2_driver_apply_fixed_step (d, &tau, A, 1, y);	
			}
			if (E >= 0)													//When energy becomes positive (i.e. escape velocity is achieved) print out the payload and terminate loop
			{
				
				fprintf(fp, "%.15f %22.15f %22.15f %22.15f\n", tb, tau, E, 1 - tau);
				printf( "%.15f %19.15f %19.15f %21.15f\n", tb, tau, E, 1 - tau);	
				break;
			}	

		}
		

		
		gsl_odeiv2_driver_free (d);
		
		tb = tb + 100;
	}
	

return 0;
}
