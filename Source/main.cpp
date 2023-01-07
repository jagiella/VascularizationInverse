/*
 * main.cpp
 *
 *  Created on: 25.01.2011
 *      Author: jagiella
 */

#include <math.h>
#include <float.h>

#include <stdio.h>
#include <stdlib.h>
#include <sys/dir.h>
#include <sys/stat.h>
#include <png.h>

#include <unistd.h>
#include <ctype.h>

#include <time.h>

//#include <conio.h>
//#include <windows.h>

#include "PNGIO.hpp"
#include "Mathematix.h"


#define MIN(a,b) (a<b?a:b)
//#define MIN(a,b,c,d) MIN(a,MIN(b,MIN(c,d)))
#define MAX(a,b) (a>b?a:b)
#define MINMAX(a,b,c) MAX(a, MIN( b, c))


#include <X11/Xlib.h>
#include <X11/keysymdef.h>

#define LOG 0
#define LIN 1
#define PARAMETER_SCALING LIN

double getParameter( double *parameters, int i)
{
#if PARAMETER_SCALING == LOG
	//return pow( 10, parameters[i]);
	return exp( parameters[i]);
#else
	return parameters[i];
#endif
}

void setParameter( double *parameters, int i, double &value)
{
#if PARAMETER_SCALING == LOG
	//parameters[i] = log10( value);
	parameters[i] = log( value);
#else
	parameters[i] = value;
#endif

}





int XGetKeyState(Display *display, KeySym keysym) {
	KeyCode keycode = XKeysymToKeycode(display, keysym);
	//Convert our key symbol to a keycode

	if (keycode == NoSymbol) {
		return 0;
	}
	//If we have an invalid keycode, abort.

	char keys[32];
	//Create a temporary array for the bits to be stored
	// possibly put this in the global scope to reduce the number
	// of times that we create it

	XQueryKeymap(display, keys);
	//Update info from the keyboard
	//If we do any long checks (multiple keys)
	//We should do this once, before the whole Loop

	return (keys[keycode >> 3] >> (keycode & 0x07)) & 0x01;
	//return the state of the choosen key
}

bool keyPressed()
{
	Display *display = XOpenDisplay(0);
	if( display==0){
		fprintf( stderr, "ERROR: Couldn't open display\n");
		exit(0);
	}else{
	//	fprintf( stderr, "Connected to display: %s\n", XDisplayString(display));
	}

	char keys_return[32];
	bool key_pressed = false;
	XQueryKeymap(display, keys_return);
	for( int i=0; i<32; i++){
		if(i!=5 && keys_return[i] > 0){
			key_pressed = true;
			/*char key = keys_return[i];
			int ii=0;
			while(key!=1){
				key=key>>1;
				ii++;
			}
			fprintf(stderr, "Key Pressed: %i, %i -> %i\n", i, ii, i*8 + ii);

			unsigned keycode = i*8 + ii;
			KeySym keysym = XKeycodeToKeysym(display, keycode, 0);
			//if(keysym == XK_space)
			//	fprintf(stderr, "Space!\n");
			//fprintf(stderr, "KeySym: %i\n", keysym);
			exit(0);
			 */
		}
	}

	XCloseDisplay(display);

	return key_pressed;
}



void Linear( int size, double *x, double *y, int degreesOfLiberty, double *parameters)
{
	double x_0 = parameters[0],
		   a   = parameters[1];
	for( int i=0; i<size; i++){
		//y[i] = a*(x[i]-x_0);
		y[i] = a*x[i]+x_0;
	}
}

void Quadratic( int size, double *x, double *y, int degreesOfLiberty, double *parameters)
{
	double x_0 = parameters[0],
		   y_0 = parameters[1],
		   a   = parameters[2];
	for( int i=0; i<size; i++){
		y[i] = a*pow(x[i]-x_0, 2.) + y_0;
	}
}

void Cubic( int size, double *x, double *y, int degreesOfLiberty, double *parameters)
{
	double a0 = parameters[0],
		   a1 = parameters[1],
		   a2 = parameters[2],
		   x_0 = parameters[3],
		   x_1 = parameters[4],
		   x_2 = parameters[5];
	for( int i=0; i<size; i++){
		y[i] = a0*pow(x[i]-x_0, 3.) + a1*pow(x[i]-x_1, 2.) + a2*pow(x[i]-x_2, 2.);
	}
}

//#define PI 3.14159265
double parkerAIF( double t,double A,double T,double sigma, double alpha,double beta,double tau,double s)
{
	return A/(sigma*sqrt(2.*PI)) * exp(-(pow(t-T,2)) / 2./pow(sigma,2.)) + alpha * exp(-beta*t) / (1.+exp(-s*(t-tau)));
}
double ArterialInputFunction( //mM
		double time//min
		)
{
	double A1 = 0.809;
	double A2 = 0.330;
	double T1 = 0.17046;
	double T2 = 0.365;
	double sigma1 = 0.0563;
	double sigma2 = 0.132;
	double alpha = 1.050;
	double beta = 0.1685;
	double s = 38.078;
	double tau = 0.483;
	double t = time;
	if( t < 0)
		return 0.;
	else{
		if(true){
			return parkerAIF(t,A1,T1,sigma1, alpha,beta,tau,s)
					+ parkerAIF(t,A2,T2,sigma2, alpha,beta,tau,s);
		}else{
			double inj=200./60.; // min
			return exp(-pow( (t-2.7*inj/4.) / (inj/4.), 2.));
		}
	}
}

double Brix2AIF( double time)
{
	return exp(-time)*(1.-exp(-time))*4.*10.;
}

void ToftsExtendedGeneralizedKineticModel_3Parameters( int size, double *t, double *y, int degreesOfLiberty, double *parameters)
{
	// Fit Parameters
	double v_P  = MINMAX( 0., parameters[0], 1.),
		   v_I  = 1. - v_P,
		   K_trans = MAX( 0., parameters[1]),
		   time0= parameters[2];

	double V = 60*60*60;
	double V_I = V*v_I;
	double V_P = V*v_P;

	// Parameters
	double dt = 0.5;

	// Variables
	double C_P = 0, C_T = 0, integral = 0;

	for( int i=0; i<size; i++){
		// integral
		integral=0;
		if( v_I>0)
		for( double tau=0; tau<=t[i]; tau+=dt){
			C_P = ArterialInputFunction( (tau-time0)/60.);
			integral+= C_P * exp( -K_trans/v_I * (t[i] - tau )) * dt;
		}

		// Update
		C_P = ArterialInputFunction( (t[i]-time0)/60.);
		y[i] = v_P*C_P + K_trans*integral;
	}
}

void Brix2( int size, double *t, double *y, int degreesOfLiberty, double *parameters)
{
	//fprintf( stderr, "Brix2()\n");
	double V_P  = MAX(0,parameters[0]),
		   V_I  = MAX(0,parameters[1]),
		   K_PS = parameters[2],
		   F    = parameters[3],
		   C_A = 1.;

	// CHECK PARAMETER INTERVALS
	if( V_P<0 || V_I<0){
		fprintf(stderr, "ERROR: Volume can not be negative: V_P=%lf, V_I=%lf", V_P, V_I);
		exit(0);
	}
	if( K_PS<0){
		fprintf(stderr, "ERROR: Exchange rate can not be negative: K_PS=%lf", K_PS);
		exit(0);
	}
	if( F<0){
		fprintf(stderr, "WARNING: Flow should not be negative: F=%lf", F);
	}

	double V = V_P + V_I;
	double v_I = (V==0?0:V_I / V);
	double v_P = (V==0?0:V_P / V);
	double dt = 0.005;
	/*if( V_I>0)
		dt = MIN(dt, MIN(fabs(K_PS/V_I), fabs(F/V_I)));
	if( V_P>0)
		dt = MIN(dt, MIN(fabs(K_PS/V_P), fabs(F/V_P)));
	fprintf( stderr, "dt=%e\n", dt);
*/
	//double dt = MIN(K_PS/V_P, F/V_P);
	//double dt = MIN(K_PS/V, F/V);

	double **A = newDoubleMatrix(2,2);
	double b[2];
	double x[2];

	double C_I = 0.,
		   C_P = 0.;

	y[0] = C_I*v_I + C_P*v_P;
	for( int i=0; i<size-1; i++){
		double time = t[i];
		y[i+1] = 0.;
		do{
			// Arterial Input Function
			//C_A = Brix2AIF( time);
			C_A = ArterialInputFunction( time/60.);
			//fprintf( stderr, "%lf\n", (F/V_P * (C_A - C_P) - K_PS/V_P * (C_P - C_I)));


			// Set linear system
			// plasma
//			C_P - C_P_old = dt * (F * (C_A - C_P) - K_PS * (C_P - C_I)) / V_P;
//			(V_P/dt + F + K_PS)*C_P - K_PS*C_I = F*C_A + C_P_old*V_P/dt;

			//dC_P = dt * (F * (C_A - C_P) - K_PS * (C_P - C_I)) / V_P;
			//C_P - C_P_old = dt * (F * (C_A - C_P) - K_PS * (C_P - C_I)) / V_P;
			//C_P/dt*V_P - C_P_old*V_P/dt = F * (C_A - C_P) - K_PS * (C_P - C_I);
			//C_P/dt*V_P - C_P_old*V_P/dt = F*C_A - F*C_P - K_PS*C_P + K_PS*C_I;
			//C_P/dt*V_P + F*C_P + K_PS*C_P - K_PS*C_I = F*C_A + C_P_old*V_P/dt;
			//C_P * ( V_P/dt + F + K_PS) - C_I * K_PS = F*C_A + C_P_old*V_P/dt;


			A[0][0] = (V_P/dt + F + K_PS);
			A[0][1] = - K_PS;
			b[0] = F*C_A + C_P*V_P/dt;

			// EES
//			dC_I = dt * (                  K_PS * (C_P - C_I)) / V_I;
//			-K_PS*C_P + (1 + K_PS)*C_I = C_I_old*V_I/dt

			//dC_I = dt * (                  K_PS * (C_P - C_I)) / V_I;
			//C_I - C_I_old = dt * (                  K_PS * (C_P - C_I)) / V_I;
			//C_I*V_I/dt - C_I_old*V_I/dt = K_PS * (C_P - C_I);
			//C_I*V_I/dt - K_PS * (C_P - C_I) = C_I_old*V_I/dt;
			//C_I * (V_I/dt + K_PS) - K_PS * C_P = C_I_old*V_I/dt;
			//
			A[1][0] = -K_PS;
			A[1][1] = (V_I/dt + K_PS);
			b[1] = C_I*V_I/dt;

			/*double scale = A[1][0]/A[0][0];
			A[1][0] = -scale*A[0][0];
			A[1][1] = -scale*A[0][1];
			b[1]    = -scale*b[0];*/
			solveLinearSystem( A, b, x, 2);

			C_P = x[0];
			C_I = x[1];



			// Plasma Concentration
			/*double dC_P = 0.;
			if(v_P>0.)
				dC_P = dt * (F * (C_A - C_P) - K_PS * (C_P - C_I)) / V_P;

			// Interstitial Space Concentration
			double dC_I = 0.;
			if(v_I>0.)
				dC_I = dt * (                  K_PS * (C_P - C_I)) / V_I;


			if(isnan(dC_P) || isnan(dC_I) || isinf(dC_P) || isinf(dC_I)){
				fprintf( stderr, "C_P=%lf C_I=%lf, C_A=%lf, t=%lf, dt=%lf\n", C_P, C_I, C_A, time, dt);
				fprintf( stderr, "V_P=%lf V_I=%lf, K_PS=%lf, F=%lf\n", V_P, V_I, K_PS, F);
				fprintf( stderr, "NUMERICAL ERROR: You should probably chose a small step size!!!\n");
				exit(0);
			}
			else{
				C_P += dC_P;
				C_I += dC_I;
			}*/
			time += dt;
		}while( time <= t[i+1]);
		y[i+1] = C_I*v_I + C_P*v_P;
	}

	deleteDoubleMatrix(A,2);
	//fprintf(stderr, " v_P=%lf v_I=%lf V_P=%lf V_I=%lf\n",v_P,v_I,V_P,V_I);
}

void Brix2_3Parameters( int size, double *t, double *y, int degreesOfLiberty, double *parameters)
{
	//fprintf( stderr, "Brix2()\n");
	/*double v_P  = MINMAX( 0., parameters[0], 1.),
		   v_I  = 1. - v_P,
		   parameters[1]),
		   F    = MAX( 0., parameters[2]),
		   C_A = 1.;*/
	double v_P  = getParameter( parameters, 0),
		   v_I  = 1. - v_P,
		   K_PS = (getParameter( parameters, 1)),
		   F    = (getParameter( parameters, 2)),
		   C_A = 1.;

	double V = 60*60*60;
	double V_I = V*v_I;
	double V_P = V*v_P;

	F=(F);

	// CHECK PARAMETER INTERVALS
	if( V_P<0 || V_I<0){
		//fprintf(stderr, "ERROR: Volume can not be negative: V_P=%lf, V_I=%lf", V_P, V_I);
		//exit(0);
		v_P  = MAX( 0., v_P);
		v_I  = MAX( 0., v_I);
	}
	if( K_PS<0){
		//fprintf(stderr, "ERROR: Exchange rate can not be negative: K_PS=%lf", K_PS);
		//exit(0);
		K_PS = MAX( 0., K_PS);
	}
	if( F<0){
		//fprintf(stderr, "WARNING: Flow should not be negative: F=%lf", F);
		F = MAX( 0., F);
	}

	double dt = 0.05; // 0.05
	/*if( V_I>0)
		dt = MIN(dt, MIN(fabs(K_PS/V_I), fabs(F/V_I)));
	if( V_P>0)
		dt = MIN(dt, MIN(fabs(K_PS/V_P), fabs(F/V_P)));
	fprintf( stderr, "dt=%e\n", dt);
*/
	//double dt = MIN(K_PS/V_P, F/V_P);
	//double dt = MIN(K_PS/V, F/V);

	double **A = newDoubleMatrix(2,2);
	double b[2];
	double x[2];

	double C_I = 0.,
		   C_P = 0.;

	y[0] = C_I*v_I + C_P*v_P;
	for( int i=0; i<size-1; i++){
		//fprintf( stderr, "t[%i] = %lf\n", i, t[i]);
		double time = t[i];
		y[i+1] = 0.;
		do{
			// Arterial Input Function
			//C_A = Brix2AIF( time);
			C_A = ArterialInputFunction( time/60. /*+ 2/60.*/);
			//fprintf( stderr, "%lf\n", (F/V_P * (C_A - C_P) - K_PS/V_P * (C_P - C_I)));


			// Set linear system
			// plasma
//			C_P - C_P_old = dt * (F * (C_A - C_P) - K_PS * (C_P - C_I)) / V_P;
//			(V_P/dt + F + K_PS)*C_P - K_PS*C_I = F*C_A + C_P_old*V_P/dt;

			//dC_P = dt * (F * (C_A - C_P) - K_PS * (C_P - C_I)) / V_P;
			//C_P - C_P_old = dt * (F * (C_A - C_P) - K_PS * (C_P - C_I)) / V_P;
			//C_P/dt*V_P - C_P_old*V_P/dt = F * (C_A - C_P) - K_PS * (C_P - C_I);
			//C_P/dt*V_P - C_P_old*V_P/dt = F*C_A - F*C_P - K_PS*C_P + K_PS*C_I;
			//C_P/dt*V_P + F*C_P + K_PS*C_P - K_PS*C_I = F*C_A + C_P_old*V_P/dt;
			//C_P * ( V_P/dt + F + K_PS) - C_I * K_PS = F*C_A + C_P_old*V_P/dt;


			A[0][0] = (V_P/dt + F + K_PS);
			A[0][1] = - K_PS;
			b[0] = F*C_A + C_P*V_P/dt;

			// EES
//			dC_I = dt * (                  K_PS * (C_P - C_I)) / V_I;
//			-K_PS*C_P + (1 + K_PS)*C_I = C_I_old*V_I/dt

			//dC_I = dt * (                  K_PS * (C_P - C_I)) / V_I;
			//C_I - C_I_old = dt * (                  K_PS * (C_P - C_I)) / V_I;
			//C_I*V_I/dt - C_I_old*V_I/dt = K_PS * (C_P - C_I);
			//C_I*V_I/dt - K_PS * (C_P - C_I) = C_I_old*V_I/dt;
			//C_I * (V_I/dt + K_PS) - K_PS * C_P = C_I_old*V_I/dt;
			//
			A[1][0] = -K_PS;
			A[1][1] = (V_I/dt + K_PS);
			b[1] = C_I*V_I/dt;

			/*double scale = A[1][0]/A[0][0];
			A[1][0] = -scale*A[0][0];
			A[1][1] = -scale*A[0][1];
			b[1]    = -scale*b[0];*/
			solveLinearSystem( A, b, x, 2);

			C_P = x[0];
			C_I = x[1];



			// Plasma Concentration
			/*double dC_P = 0.;
			if(v_P>0.)
				dC_P = dt * (F * (C_A - C_P) - K_PS * (C_P - C_I)) / V_P;

			// Interstitial Space Concentration
			double dC_I = 0.;
			if(v_I>0.)
				dC_I = dt * (                  K_PS * (C_P - C_I)) / V_I;


			if(isnan(dC_P) || isnan(dC_I) || isinf(dC_P) || isinf(dC_I)){
				fprintf( stderr, "C_P=%lf C_I=%lf, C_A=%lf, t=%lf, dt=%lf\n", C_P, C_I, C_A, time, dt);
				fprintf( stderr, "V_P=%lf V_I=%lf, K_PS=%lf, F=%lf\n", V_P, V_I, K_PS, F);
				fprintf( stderr, "NUMERICAL ERROR: You should probably chose a small step size!!!\n");
				exit(0);
			}
			else{
				C_P += dC_P;
				C_I += dC_I;
			}*/
			time += dt;
		}while( time <= t[i+1]);
		y[i+1] = C_I*v_I + C_P*v_P;
	}

	deleteDoubleMatrix(A,2);
	//fprintf(stderr, " v_P=%lf v_I=%lf V_P=%lf V_I=%lf\n",v_P,v_I,V_P,V_I);
}

void Brix2_4Parameters_Delay( int size, double *t, double *y, int degreesOfLiberty, double *parameters)
{
	//fprintf( stderr, "Brix2()\n");
	double v_P  = MINMAX( 0., parameters[0], 1.),
		   v_I  = 1. - v_P,
		   K_PS = MAX( 0., parameters[1]),
		   F    = MAX( 0., parameters[2]),
		   time0= parameters[3],
		   C_A = 1.;

	double V = 60*60*60;
	double V_I = V*v_I;
	double V_P = V*v_P;

	F=(F);

	//K_PS = V_P * K_PS;

	// CHECK PARAMETER INTERVALS
	if( V_P<0 || V_I<0){
		fprintf(stderr, "ERROR: Volume can not be negative: V_P=%lf, V_I=%lf", V_P, V_I);
		exit(0);
	}
	if( K_PS<0){
		fprintf(stderr, "ERROR: Exchange rate can not be negative: K_PS=%lf", K_PS);
		exit(0);
	}
	if( F<0){
		fprintf(stderr, "WARNING: Flow should not be negative: F=%lf", F);
	}

	double dt = 0.5;
	/*if( V_I>0)
		dt = MIN(dt, MIN(fabs(K_PS/V_I), fabs(F/V_I)));
	if( V_P>0)
		dt = MIN(dt, MIN(fabs(K_PS/V_P), fabs(F/V_P)));
	fprintf( stderr, "dt=%e\n", dt);
*/
	//double dt = MIN(K_PS/V_P, F/V_P);
	//double dt = MIN(K_PS/V, F/V);

	double **A = newDoubleMatrix(2,2);
	double b[2];
	double x[2];

	double C_I = 0.,
		   C_P = 0.;

	//int steps = (int) ((t[size-1] - t[0])/dt);
	//double *C_P_history = (double*) malloc(steps * sizeof(double));

	y[0] = C_I*v_I + C_P*v_P;
	for( int i=0; i<size-1; i++){
		double time = t[i];
		y[i+1] = 0.;
		do{
			// Arterial Input Function
			//C_A = Brix2AIF( time);
			//if(i==0)
				C_A = ArterialInputFunction( (time-time0)/60.);
			//else
			//	C_A = C_P_history[(int)(time/dt)];

			//fprintf( stderr, "%lf\n", (F/V_P * (C_A - C_P) - K_PS/V_P * (C_P - C_I)));


			// Set linear system
			// plasma
//			C_P - C_P_old = dt * (F * (C_A - C_P) - K_PS * (C_P - C_I)) / V_P;
//			(V_P/dt + F + K_PS)*C_P - K_PS*C_I = F*C_A + C_P_old*V_P/dt;

			//dC_P = dt * (F * (C_A - C_P) - K_PS * (C_P - C_I)) / V_P;
			//C_P - C_P_old = dt * (F * (C_A - C_P) - K_PS * (C_P - C_I)) / V_P;
			//C_P/dt*V_P - C_P_old*V_P/dt = F * (C_A - C_P) - K_PS * (C_P - C_I);
			//C_P/dt*V_P - C_P_old*V_P/dt = F*C_A - F*C_P - K_PS*C_P + K_PS*C_I;
			//C_P/dt*V_P + F*C_P + K_PS*C_P - K_PS*C_I = F*C_A + C_P_old*V_P/dt;
			//C_P * ( V_P/dt + F + K_PS) - C_I * K_PS = F*C_A + C_P_old*V_P/dt;


			A[0][0] = (V_P/dt + F + K_PS);
			A[0][1] = - K_PS;
			b[0] = F*C_A + C_P*V_P/dt;

			// EES
//			dC_I = dt * (                  K_PS * (C_P - C_I)) / V_I;
//			-K_PS*C_P + (1 + K_PS)*C_I = C_I_old*V_I/dt

			//dC_I = dt * (                  K_PS * (C_P - C_I)) / V_I;
			//C_I - C_I_old = dt * (                  K_PS * (C_P - C_I)) / V_I;
			//C_I*V_I/dt - C_I_old*V_I/dt = K_PS * (C_P - C_I);
			//C_I*V_I/dt - K_PS * (C_P - C_I) = C_I_old*V_I/dt;
			//C_I * (V_I/dt + K_PS) - K_PS * C_P = C_I_old*V_I/dt;
			//
			A[1][0] = -K_PS;
			A[1][1] = (V_I/dt + K_PS);
			b[1] = C_I*V_I/dt;

			/*double scale = A[1][0]/A[0][0];
			A[1][0] = -scale*A[0][0];
			A[1][1] = -scale*A[0][1];
			b[1]    = -scale*b[0];*/
			solveLinearSystem( A, b, x, 2);

			if( V_P != 0){
				C_P = x[0];
				C_I = x[1];
			}



			// Plasma Concentration
			/*double dC_P = 0.;
			if(v_P>0.)
				dC_P = dt * (F * (C_A - C_P) - K_PS * (C_P - C_I)) / V_P;

			// Interstitial Space Concentration
			double dC_I = 0.;
			if(v_I>0.)
				dC_I = dt * (                  K_PS * (C_P - C_I)) / V_I;


			if(isnan(dC_P) || isnan(dC_I) || isinf(dC_P) || isinf(dC_I)){
				fprintf( stderr, "C_P=%lf C_I=%lf, C_A=%lf, t=%lf, dt=%lf\n", C_P, C_I, C_A, time, dt);
				fprintf( stderr, "V_P=%lf V_I=%lf, K_PS=%lf, F=%lf\n", V_P, V_I, K_PS, F);
				fprintf( stderr, "NUMERICAL ERROR: You should probably chose a small step size!!!\n");
				exit(0);
			}
			else{
				C_P += dC_P;
				C_I += dC_I;
			}*/

			//C_P_history[(int)(time/dt)] = C_P;
			time += dt;
		}while( time <= t[i+1]);
		y[i+1] = C_I*v_I + C_P*v_P;
	}

	//free(C_P_history);
	deleteDoubleMatrix(A,2);
	//fprintf(stderr, " v_P=%lf v_I=%lf V_P=%lf V_I=%lf\n",v_P,v_I,V_P,V_I);
}

void Brix2NoPermeability_2Parameters_Delay( int size, double *t, double *y, int degreesOfLiberty, double *parameters)
{
	//fprintf( stderr, "Brix2()\n");
	double v_P  = MINMAX( 0., parameters[0], 1.),
		    time0= MAX( 0., parameters[1]);

	for( int i=0; i<size; i++){
		double C_A = ArterialInputFunction( (t[i]-time0)/60.);
		y[i] = v_P*C_A;
	}
}

/*void Brix2NoPermeability_2Parameters_Flow( int size, double *t, double *y, int degreesOfLiberty, double *parameters)
{
	//fprintf( stderr, "Brix2()\n");
	double v_P = MINMAX( 0., parameters[0], 1.),
		    F   = MAX( 0., parameters[1]);

	for( int i=0; i<size; i++){
		double C_A = ArterialInputFunction( (t[i]-v_P*V/F)/60.);
		y[i] = v_P*C_A;
	}
}*/

double SumOfSquares( int size, double *data1, double *data2)
{
	double n = 2.;
	double S = 0.;
	for( int i=0; i<size; i++)
		S += pow(data1[i] - data2[i], n);
	return pow( S, 1./n);
	return S;
}

void getGradient(
		int size, double *sampleX, double *sampleY,
		int parameterSize, double *parameters, double *gradient,
		void (*f)(int, double*, double*, int, double *),
		double dx)
{
	double tempY[size];

	for( int i=0; i<parameterSize; i++)
	{
		// the change of least squares in respect to

		// RIGHT
		/*parameters[i] += dx;
		f( size, t, data2, degreesOfLiberty, parameters);
		gradient[i] = (S - leastSquareFit( size, data1, data2))/dx;
		parameters[i] -= dx;*/

		// RIGHT + LEFT
		parameters[i] += dx;

		f( size, sampleX, tempY, parameterSize, parameters);
		gradient[i] = SumOfSquares( size, sampleY, tempY);

		parameters[i] -= 2*dx;

		f( size, sampleX, tempY, parameterSize, parameters);
		gradient[i] -=SumOfSquares( size, sampleY, tempY);
		gradient[i] /= 2.*dx;

		parameters[i] += dx;


		// LEFT
		/*parameters[i] -= dx;
		f( size, t, data2, degreesOfLiberty, parameters);
		gradient[i] = (leastSquareFit( size, data2, data1) - S)/dx;
		parameters[i] += dx;
		 */
	}
}

void getJacobiMatrix(
		int size, double *sampleX, double *sampleY,
		int parameterSize, double *parameters, double **jacobi,
		void (*f)(int, double*, double*, int, double *),
		double dx)
{
	double tempY[size];

	for( int j=0; j<parameterSize; j++)
	{
		// the change of least squares in respect to

		// RIGHT
		/*parameters[i] += dx;
		f( size, t, data2, degreesOfLiberty, parameters);
		gradient[i] = (S - leastSquareFit( size, data1, data2))/dx;
		parameters[i] -= dx;*/

		// RIGHT + LEFT
		parameters[j] += dx;

		//fprintf(stderr, "before\n");
		f( size, sampleX, tempY, parameterSize, parameters);
		//fprintf(stderr, "after\n");

		for( int i=0; i<size; i++){
			jacobi[i][j] = tempY[i];
		}
		parameters[j] -= dx;

		parameters[j] -= dx;

		//fprintf(stderr, "before2\n");
		f( size, sampleX, tempY, parameterSize, parameters);
		//fprintf(stderr, "after2\n");
		for( int i=0; i<size; i++){
			jacobi[i][j] -= tempY[i];
			jacobi[i][j] /= 2.*dx;
		}

		parameters[j] += dx;


		// LEFT
		/*parameters[i] -= dx;
		f( size, t, data2, degreesOfLiberty, parameters);
		gradient[i] = (leastSquareFit( size, data2, data1) - S)/dx;
		parameters[i] += dx;
		 */
	}
}

double randAB( double a, double b)
{
	return a + (b-a) * rand() / RAND_MAX;
}

double randNormalAB( double mean, double std)
{
	double u1, u2, q, p;	
	do{
		u1 = randAB( -1, +1);
		u2 = randAB( -1, +1);

		q = u1*u1 + u2*u2;
	}while( q == 0 || q > 1);

	p = sqrt(-2*log(q)/q);

	return mean + std*p*u1;
}

void LevenbergMarquardt(
	int sampleSize, double *sampleX, double *sampleY,
	int parameterSize, double *parameters, double *parametersMin, double *parametersMax,
	double differentiationStepSize, double lambda, double minSquares, int maxIterations,
	void (*f)(int, double*, double*, int, double*)
	)
{
	fprintf(stderr, "Start LevenbergMarquardt\n");

	double S, S_old = 0.;
	int it=0;
	double **jacobi = newDoubleMatrix(sampleSize, parameterSize);
	double fitY[sampleSize];

	double **A = newDoubleMatrix(parameterSize, parameterSize);
	double b[parameterSize];
	double x[parameterSize];

	// COPY PARAMETERS
	double log_parameters[parameterSize];
	double log_parametersMin[parameterSize];
	double log_parametersMax[parameterSize];
	for (int i = 0; i < parameterSize; i++) {
		setParameter( log_parameters,    i, parameters[i]);
		setParameter( log_parametersMin, i, parametersMin[i]);
		setParameter( log_parametersMax, i, parametersMax[i]);
	}


	// get function sample
	(*f)( sampleSize, sampleX, fitY, parameterSize, log_parameters);

	// compare data to fit
	S = SumOfSquares( sampleSize, sampleY, fitY);
	fprintf(stderr, "0 iterations: S = %e -> step size = %e \b", S, differentiationStepSize);



	for( ; /*!keyPressed() &&*/ it<maxIterations /*&& fabs(S-S_old) > 1e-15*/ && S>minSquares; it++)
	{

		// Construct Jacobi-Matrix
		getJacobiMatrix(sampleSize, sampleX, sampleY, parameterSize, log_parameters, jacobi, f, differentiationStepSize);

		// Construct Linear System
		for (int i = 0; i < parameterSize; i++) {

			// Levenberg
			for (int j = 0; j < parameterSize; j++) {
				A[i][j] = 0.;
				for (int k = 0; k < sampleSize; k++)
					// A = J^T J
					A[i][j] += jacobi[k][i] * jacobi[k][j];
			}
			// Marquardt
			A[i][i] *= (1. + lambda);

			b[i] = 0.;
			for (int k = 0; k < sampleSize; k++)
				b[i] -= jacobi[k][i] * (sampleY[k] - fitY[k]);
		}

		// Solve Linear System
		solveLinearSystem(A, b, x, parameterSize);

		// Update Parameters
		for (int i = 0; i < parameterSize; i++)
			log_parameters[i] -= x[i];

		// Check upper and lower bounds
		for( int i=0; i<parameterSize; i++){
			if( log_parameters[i] < log_parametersMin[i]) log_parameters[i] = log_parametersMin[i];
			if( log_parameters[i] > log_parametersMax[i]){
				log_parameters[i] = log_parametersMax[i];
				//fprintf( stderr, "UPPER BOUND: log(%e) = %e\n", exp(log_parameters[i]), log_parameters[i]);
			}
		}


		// get function sample
		(*f)( sampleSize, sampleX, fitY, parameterSize, log_parameters);

		// squares
		S_old = S;
		S = SumOfSquares( sampleSize, sampleY, fitY);
		fprintf(stderr, "\r%i iterations: S = %e -> step size = %e \b", it, S, differentiationStepSize);
	}

	// COPY BACK PARAMETERS
	for (int i = 0; i < parameterSize; i++) {
		parameters[i]    = getParameter( log_parameters, i);
		parametersMin[i] = getParameter( log_parametersMin, i);
		parametersMax[i] = getParameter( log_parametersMax, i);
	}



	deleteDoubleMatrix( A, parameterSize);
	deleteDoubleMatrix( jacobi, sampleSize);

	fprintf(stderr, "\nParameters: ");
	for( int i=0; i<parameterSize; i++)
		fprintf(stderr, "%10.3e ", parameters[i]);
	fprintf(stderr, "\n");
}


void GradientDecent(
		int sampleSize, double *sampleX, double *sampleY,
		int parameterSize, double *parameters, double *parametersMin, double *parametersMax,
		double differentiationStepSize, double gamma, double minSquares, int maxIterations,
		void (*f)(int, double*, double*, int, double*)
		)
	{
	// FIT
	double S, S_old = 0.;
	int it=0;
	int decent = 0;
	double fitY[sampleSize];
	double gradient[parameterSize];


	//

	// get function sample
	(*f)( sampleSize, sampleX, fitY, parameterSize, parameters);

	// compare data to fit
	S = SumOfSquares( sampleSize, sampleY, fitY);
	fprintf(stderr, "0 iterations: S = %lf -> gamma = %lf\n", S, gamma);

	for( ; !keyPressed() && it<maxIterations && fabs(S-S_old) > 1e-10 && S>minSquares; it++)
	{

		// get gradient of least squares
		//fprintf(stderr, "gradient\n");
		double differentiationStepSize = 0.0001;
		getGradient(
				sampleSize, sampleX, sampleY,
				parameterSize, parameters, gradient,
				f,
				differentiationStepSize);
		/*fprintf(stderr, "Parameters: ");
		for( int i=0; i<parameterSize; i++)
			fprintf(stderr, "%10.3lf ", parameters[i]);
		fprintf(stderr, "\n");
		fprintf(stderr, "Gradient:   ");
		for( int i=0; i<parameterSize; i++)
			fprintf(stderr, "%10.3lf ", gradient[i]);
		fprintf(stderr, "\n");*/
		//fprintf(stderr, "end gradient\n");

		// update parameters
		for( int i=0; i<parameterSize; i++){
			/*if(parameters[i] - stepSize*gradient[i] < parametersMin[i]){
				fprintf( stderr, "Parameter[%i] too small: %lf\n", i, parameters[i] - stepSize*gradient[i]);
				gradient[i] = (parameters[i] - parametersMin[i])/stepSize;
			}
			if(parameters[i] - stepSize*gradient[i] > parametersMax[i]){
				fprintf( stderr, "Parameter[%i] too large: %lf\n", i, parameters[i] - stepSize*gradient[i]);
				gradient[i] = (parameters[i] - parametersMax[i])/stepSize;
			}*/
			if( parameters[i] - gamma*gradient[i] < parametersMin[i]){
				gradient[i] = (parameters[i] - parametersMin[i])/gamma;
				parameters[i] = parametersMin[i];
			}else if( parameters[i] - gamma*gradient[i] > parametersMax[i]){
				gradient[i] = (parameters[i] - parametersMax[i])/gamma;
				parameters[i] = parametersMax[i];
			}else
				parameters[i] -= gamma*gradient[i];
			/*if(parametersMin != 0){
				parameters[i] = MAX(parameters[i], parametersMin[i]);
			}
			if(parametersMax != 0)
				parameters[i] = MIN(parameters[i], parametersMax[i]);*/
		}


		// get function sample
		(*f)( sampleSize, sampleX, fitY,
			  parameterSize, parameters);

		// squares
		S_old = S;
		S = SumOfSquares( sampleSize, sampleY, fitY);
		fprintf(stderr, "%i iterations: S = %lf -> gamma = %lf\n", it, S, gamma);

		// change acception & step size adaptation
		if(S > S_old){
			// reverse parameter change
			for( int i=0; i<parameterSize; i++)
				parameters[i] += gamma*gradient[i];


			decent = 0;
			gamma /= 2.;

		}
		else{
			decent ++;
			if(decent >= 5)
				gamma *= 1.1;
		}

		/*for( int i=0; i<parameterSize; i++){
	//		parameters[i] -= stepSize*gradient[i];
			if(parametersMin != 0)
				parameters[i] = MAX(parameters[i], parametersMin[i]);
			if(parametersMax != 0)
				parameters[i] = MIN(parameters[i], parametersMax[i]);
		}*/

	}

	(*f)( sampleSize, sampleX, fitY,
		  parameterSize, parameters);

	fprintf(stderr, "Parameters: ");
	for( int i=0; i<parameterSize; i++)
		fprintf(stderr, "%10.3lf ", parameters[i]);
	fprintf(stderr, "\n");
	S = SumOfSquares( sampleSize, sampleY, fitY);

	/*fprintf(stderr, "Fit after %i iterations (error=%e)\n", it, S);
	for( int i=0; i<sampleSize; i++){
		fprintf(stderr, "%10lf -> %10lf\n", sampleX[i], fitY[i]);
	}*/
}

/*
int color_type;
int bit_depth;
png_structp png_ptr;
png_infop info_ptr;
png_uint_32 width, height;

int readpng_init(FILE *infile, long *pWidth, long *pHeight)
{
	unsigned char sig[8];

	fread(sig, 1, 8, infile);
	if (!png_check_sig(sig, 8))
		return 1;   // bad signature

	png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	if (!png_ptr)
		return 4;   // out of memory

	info_ptr = png_create_info_struct(png_ptr);
	if (!info_ptr) {
		png_destroy_read_struct(&png_ptr, NULL, NULL);
		return 4;   // out of memory
	}

	if (setjmp(png_ptr->jmpbuf)) {
		png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
		return 2;
	}

	png_init_io(png_ptr, infile);
	png_set_sig_bytes(png_ptr, 8);
	png_read_info(png_ptr, info_ptr);

	png_get_IHDR(png_ptr, info_ptr, &width, &height, &bit_depth, &color_type, NULL, NULL, NULL);
	*pWidth = width;
	*pHeight = height;

	return 0;
}

unsigned char *readpng_get_image(double display_exponent, int *pChannels, unsigned long *pRowbytes)
{
	if (color_type == PNG_COLOR_TYPE_PALETTE)
		png_set_expand(png_ptr);
	if (color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8)
		png_set_expand(png_ptr);
	if (png_get_valid(png_ptr, info_ptr, PNG_INFO_tRNS))
		png_set_expand(png_ptr);

	// These functions are FICTITIOUS!  They DO NOT EXIST in any
	// version of libpng to date (through 1.0.3).

	if (color_type == PNG_COLOR_TYPE_PALETTE)
		png_set_palette_to_rgb(png_ptr);
	if (color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8)
		png_set_gray_1_2_4_to_8(png_ptr);
	if (png_get_valid(png_ptr, info_ptr, PNG_INFO_tRNS))
		png_set_tRNS_to_alpha(png_ptr);

	  if (bit_depth == 16)
		png_set_strip_16(png_ptr);
	if (color_type == PNG_COLOR_TYPE_GRAY || color_type
			== PNG_COLOR_TYPE_GRAY_ALPHA)
		png_set_gray_to_rgb(png_ptr);

	double gamma;

	if (png_get_gAMA(png_ptr, info_ptr, &gamma))
		png_set_gamma(png_ptr, display_exponent, gamma);

	png_uint_32 i, rowbytes;
	png_bytep row_pointers[height];

	png_read_update_info(png_ptr, info_ptr);

	fprintf(stderr, "test1\n");
	*pRowbytes = rowbytes = png_get_rowbytes(png_ptr, info_ptr);
	fprintf(stderr, "test2\n");
	*pChannels = (int) png_get_channels(png_ptr, info_ptr);
	fprintf(stderr, "test3\n");

	unsigned char *image_data;
	if ((image_data = (unsigned char *) malloc(rowbytes * height)) == NULL) {
		png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
		return NULL;
	}

	for (i = 0; i < height; ++i)
		//row_pointers[i] = image_data + i * rowbytes;
		row_pointers[i] = &image_data[i * rowbytes];

	fprintf(stderr, "test4\n");
	png_read_image(png_ptr, row_pointers);

	fprintf(stderr, "test5\n");
	png_read_end(png_ptr, NULL);

	fprintf(stderr, "test6\n");
	return image_data;
}
*/

#define BRIX2		0
#define BRIX2DELAY	1
#define BRIX2NOPERM 2
#define TOFTS 		3
#define TOFTSNODELAY	4



int main( int argc, char **argv)
{
	// READ DATA FROM FILE
	/*FILE *fp_png = fopen("test.png","rb");
	long pWidth;
	long pHeight;
	fprintf( stderr, "png: %i\n", readpng_init(fp_png, &pWidth, &pHeight));
	fclose(fp_png);

	fprintf( stderr, "image -> %li x %li\n", pWidth, pHeight);

	double CRT_exponent = 2.2;
	double LUT_exponent = 1.8 / 2.61;
	double default_display_exponent = LUT_exponent * CRT_exponent;
	int channels = 0;
	unsigned long rowbytes = 0;
	unsigned char *image_data = readpng_get_image(default_display_exponent, &channels, &rowbytes);
*/

	/*PngIO *pngImg = new PngIO();

	pngImg->read_png_file( (char*)"test.png");
	pngImg->process_file();
	pngImg->readpng_cleanup(0);
	pngImg->write_png_file((char*)"test2.png");
	pngImg->writepng_cleanup(1);

	delete( pngImg);*/

	// OPTIONS
	char c;
	//int aflag = 0;
	//int bflag = 0;
	//char *cvalue = NULL;
	double scale_x = 1.;
	double scale_y = 1.;// 1/256.;
	int coarsing = 1;

	int pixelX=-1;
	int pixelY=-1;

	int    Imax = 1;
	double Smin = 0.1;
	double stepSize = 1;
	double fitMethodParameter = 1;
	double noise = 0;

	int inverseMethode = BRIX2;

	bool ROI = false;

	char dirname[512] = ".";

	int index;
	while ((c = getopt(argc, argv, "D:d:a:c:n:x:y:X:Y:S:I:M:R")) != -1)
		switch (c) {
		case 'n':{
			noise = atof(optarg);
		}break;
		case 'c':{
			coarsing = atoi(optarg);
		}break;
		case 'D':{

			sprintf( dirname, "%s/", optarg);
						if( mkdir(dirname, S_IRWXU)!=0){
							fprintf(stderr, "WARNING: Can not create directory %s: ", optarg);
							perror("");
							//perror("mkdir() error");
							//exit(0);
						}else{
							fprintf(stderr, "Create Directory %s\n", optarg);
						}

						FILE *fp;
						char filename[512];
						sprintf( filename, "%s/commandline.dat", dirname);
						fprintf(stderr, "write command line to %s\n", filename);
						fp = fopen( filename, "a+");
						if( !fp){
							fprintf(stderr, "Couldn't create file\n");
							perror("sd");
						}
						time_t rawtime;
						time( &rawtime );
						fprintf( fp, "\n%s", asctime(localtime(&rawtime)));
						for( int i=0; i<argc; i++)
							fprintf( fp, "%s ", argv[i]);
						fprintf( fp, "\n");
						fclose(fp);

		}break;
		case 'M':
			if( strstr(optarg, "BRIX2NOPERM") != 0){
				inverseMethode = BRIX2NOPERM;
			}else if( strstr(optarg, "BRIX2DELAY") != 0){
				inverseMethode = BRIX2DELAY;
			}else if( strstr(optarg, "BRIX2") != 0){
				inverseMethode = BRIX2;
			}else if( strstr(optarg, "TOFTSNODELAY") != 0){
				inverseMethode = TOFTSNODELAY;
			}else if( strstr(optarg, "TOFTS") != 0){
				inverseMethode = TOFTS;
			}
			fprintf(stderr, "Fitting Function is %s.\n", optarg);
			break;
		case 'd':
			stepSize = atof(optarg);
			fprintf(stderr, "Numerically solve ODE's with step size h = %e.\n", stepSize);
			break;
		case 'a':
			fitMethodParameter = atof(optarg);
			fprintf(stderr, "Fit Method specific parameter (GD: gamma, LM: lambda) = %e.\n", fitMethodParameter);
			break;
		case 'x':
			scale_x = atof(optarg);
			fprintf(stderr, "Scale in x-direction: x*%lf.\n", scale_x);
			break;
		case 'y':
			scale_y = atof(optarg);
			fprintf(stderr, "Scale in y-direction: y*%lf.\n", scale_y);
			//scale_y; //= 256.;
			break;
		case 'X':
			pixelX = atoi(optarg);
			break;
		case 'Y':
			pixelY = atoi(optarg);
			break;
		case 'I':
			Imax = atoi(optarg);
			fprintf(stderr, "Accept fits after max %i iterations.\n", Imax);
			break;
		case 'S':
			Smin = atof(optarg);
			fprintf(stderr, "Accept fits with least square values of < %e.\n", Smin);
			break;
		case 'R':
			ROI = true;
			fprintf(stderr, "Use Region of Interest (ROI: 1...i)\n");
			break;
		case '?':
			if (optopt == 'c')
				fprintf(stderr, "Option -%c requires an argument.\n", optopt);
			else if (isprint(optopt))
				fprintf(stderr, "Unknown option `-%c'.\n", optopt);
			else
				fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
			return 1;
		default:
			abort();
		}

	// INPUT
	bool png_input = false;
	bool binary_input = true;

	// TEST
	if(false){
		PngIO *pngImg = new PngIO();
		//for (index = optind; index < argc; index++){
			//printf("Non-option argument %s\n", argv[index]);

			if (pngImg->read_png_file(argv[optind])) {

					for (int x = 0; x < pngImg->getWidth(); x++){
						for (int y = 0; y < pngImg->getHeight(); y++) {
							fprintf( stderr, "(%3i,%3i,%3i) ", pngImg->getPixel(x, y, 0), pngImg->getPixel(x, y, 1), pngImg->getPixel(x, y, 2));
						}
						fprintf( stderr, "\n");
					}
			}
			pngImg->readpng_cleanup(1);
		//}
		delete( pngImg);
		exit(0);
	}
	// TEST END

	//printf("aflag = %d, bflag = %d, cvalue = %s\n", aflag, bflag, cvalue);

	// ANALYSE IMAGE FILES
	int countImageFiles = 0;
	int imageHeight = 0;
	int imageWidth = 0;
	int imageDepth = 0;
	if(png_input)
	for (index = optind; index < argc; index++){
		// for all arguments
		// is file?
		PngIO *pngImg = new PngIO();
		if( pngImg->read_png_file( argv[index])){
			countImageFiles++;

			// has good size?
			if (countImageFiles == 1) {
				imageHeight = pngImg->getHeight();
				imageWidth = pngImg->getWidth();
			} else {
				if (imageHeight != pngImg->getHeight() || imageWidth
						!= pngImg->getWidth()) {
					fprintf(
							stderr,
							"ERROR: Dimensions of image %s (%ix%i) differes from other images (%ix%i)\n",
							argv[index], pngImg->getWidth(),
							pngImg->getHeight(), imageWidth, imageHeight);
					exit(0);
				}
			}
		}else{
			fprintf(stderr, "WARNING: Can not read file %s\n", argv[index]);
		}

		pngImg->readpng_cleanup(1);
		delete( pngImg);
	}else if(binary_input){
		for (index = optind; index < argc; index++){
			FILE *fp;
			if( (fp = fopen( argv[index], "rb+")) != 0){

				countImageFiles++;

				fread( &imageWidth, sizeof(int), 1, fp);
				fread( &imageHeight, sizeof(int), 1, fp);
				fread( &imageDepth, sizeof(int), 1, fp);
				if( imageDepth>1){
					fprintf(stderr, "ATTENTION: binary data is 3-dimensional (ignoring 3rd dimension)!\n");
				}

				fclose(fp);
			}else{
				fprintf(stderr, "WARNING: Can not read file %s\n", argv[index]);
			}
		}
	}

	fprintf( stderr, "%i files found [%i x %i x %i]\n", countImageFiles,
			imageWidth,imageHeight,imageDepth);

	//exit(0);

	int sampleSize = argc - optind;
	//fprintf( stderr, "sample size = %i\n", sampleSize); exit(0);
	double sampleX[sampleSize];
	double sampleY[sampleSize];
	double ***samplesX = (double***) malloc(sizeof(double**) * (imageWidth));
	double ***samplesY = (double***) malloc(sizeof(double**) * (imageWidth));
	int    samples[imageWidth][imageHeight];
	for( int x=0; x<(imageWidth); x++){
		samplesX[x] = (double**) malloc(sizeof(double*) * imageHeight);
		samplesY[x] = (double**) malloc(sizeof(double*) * imageHeight);
		for( int y=0; y<(imageHeight); y++){
			samplesX[x][y] = (double*) malloc(sizeof(double) * countImageFiles);
			samplesY[x][y] = (double*) malloc(sizeof(double) * countImageFiles);
		}
	}

	for (int x = 0; x < imageWidth; x++)
		for (int y = 0; y < imageHeight; y++){
			samples[x][y]=0;
			for (int i = 0; i < countImageFiles; i++)
			{
				samplesX[x][y][i]=0;
				samplesY[x][y][i]=0;
			}
		}


	// NON-OPTION ARGUMENTS
	if( png_input){
		PngIO *pngImg = new PngIO();
		for (index = optind; index < argc; index++){
			//printf("Non-option argument %s\n", argv[index]);

			if (pngImg->read_png_file(argv[index])) {

				if (pixelX >= 0 && pixelY >= 0) {
					sampleX[index - optind] = scale_x * (index - optind);
					sampleY[index - optind] = scale_y * (double)pngImg->getPixelIntensity(pixelX, pixelY, 0);
				} else
					for (int x = 0; x < pngImg->getWidth(); x++)
						for (int y = 0; y < pngImg->getHeight(); y++) {
							samplesX[x][y][index - optind] = scale_x * (index - optind);
							//samplesY[x][y][index - optind] = scale_y * pngImg->getPixel(x, y, 0);
							samplesY[x][y][index - optind] = scale_y * (double)pngImg->getPixelIntensity(x, y, 0);

							//if (x == 10 && y == 10)
							//	fprintf(stderr, "(%i, %i) = %i\n", x, y, pngImg->getPixel(x, y, 0));
						}
			}
			pngImg->readpng_cleanup(1);
		}
		delete( pngImg);
	}else if( binary_input){
		FILE *fp;
		float valueold=0;
		for (index = optind; index < argc; index++){
			//printf("Non-option argument %s\n", argv[index]);

			if ((fp = fopen( argv[index], "rb+")) != 0) {

				// READ HEADER
				fread( &imageWidth, sizeof(int), 1, fp);
				fread( &imageHeight, sizeof(int), 1, fp);
				fread( &imageDepth, sizeof(int), 1, fp);

				for (int x = 0; x < imageWidth; x++)
					for (int y = 0; y < imageHeight; y++)
						for (int z = 0; z < imageDepth; z++) {

							// READ VALUE FROM FILE
							float valuef;
							fread( &valuef, sizeof(float), 1, fp);

							// ADD NOISE
							//valuef *= (1. + noise * randAB( -1, 1));
							valuef = randNormalAB( valuef, noise*valuef);

							if( ROI){
								//valueold+=valuef;
								//valuef=valueold;
								if(x>1)
								valuef = (valuef + (float)(x-1)*samplesY[x-1][y][index - optind]/scale_y)/(float)(x+1-1);
							}

							if( coarsing>1 || z==imageDepth/2) // skip 3d data
							{
								if (pixelX >= 0 && pixelY >= 0) {
									// SINGLE VOXEL
									if( x==pixelX && y==pixelY){
										sampleX[index - optind] = scale_x * (index - optind);
										sampleY[index - optind] = scale_y * valuef;
									}
								}else{
									// VOXEL MAP
									samplesX[x/coarsing][y/coarsing][index - optind] += scale_x * (index - optind);
									samplesY[x/coarsing][y/coarsing][index - optind] += scale_y * valuef;

									if(index == optind)
										samples[x/coarsing][y/coarsing]++;
								}
							}
					}

				fclose(fp);


			}
		}
	}

	// COARSING

	if( coarsing>1){
		imageWidth =  imageWidth/coarsing;
		imageHeight =  imageHeight/coarsing;

		for( int i=0; i<sampleSize; i++)
		for (int x = 0; x < imageWidth; x++) {
			for (int y = 0; y < imageHeight; y++) {
				samplesX[x][y][i] /= samples[x][y];
				samplesY[x][y][i] /= samples[x][y];
			}
		}
	}
	//END COARSING

	//return(0);




	// INIT DATA
	/*fprintf(stderr, "Data\n");
	for( int i=0; i<sampleSize; i++){
		sampleX[i] = (double)i;
		sampleY[i] = randAB(0,5);
		sampleY[i] = MAX(0, exp(-i/5.)*(1.-exp(-i/5.))*4. + 1.) * (1. + randAB(-0.05,0.05));

		//sampleY[i] = (2.*pow(i-10., 2.)+10) * (1. + randAB(-0.05,0.05));
		fprintf(stderr, "%10lf -> %10lf\n", sampleX[i], sampleY[i]);
	}*/


	// INIT PARAMETERS

	void (*f)(int, double*, double*, int, double *);
	int parameterSize;
	char parameterNames[4][512] = {"v_P",	"K_PS",		"F",		"t0"};
	double parameters[4]        = {0.0139, 	0, 		1e5, 		0};
	double parametersMin[4]     = {0, 	0,     		0,		0};
	double parametersMax[4]     = {1., 	1e10, 		1e14, 		1000};
	/*char parameterNames[2][512] = {"v_P",	"t0"};
	double parameters[2]        = {0.001, 	1};
	double parametersMin[2]     = {0., 		0};
	double parametersMax[2]     = {1., 		1000};*/

	switch( inverseMethode){
	case BRIX2DELAY:
		f = Brix2_4Parameters_Delay;
		parameterSize = 4;
		break;

	case BRIX2:
		f = Brix2_3Parameters;
		parameterSize = 3;
		/*for( int p=1; p<3; p++){
			parameters[p] = log(parameters[p]);
			parametersMin[p] = log(parametersMin[p]);
			parametersMax[p] = log(parametersMax[p]);
		}*/
		break;

	case BRIX2NOPERM:
		f = Brix2NoPermeability_2Parameters_Delay;
		parameterSize = 2;
		sprintf( parameterNames[1], "t0");
		parametersMax[1]=DBL_MAX;
		break;
	case TOFTS:
		f = ToftsExtendedGeneralizedKineticModel_3Parameters;
		parameterSize = 3;
		sprintf( parameterNames[1], "K_trans");
		sprintf( parameterNames[2], "t0");
		parametersMax[2]=DBL_MAX;
		parameters[2]=1e-6;
		break;
	case TOFTSNODELAY:
		f = ToftsExtendedGeneralizedKineticModel_3Parameters;
		parameterSize = 3;
		sprintf( parameterNames[1], "K_trans");
		sprintf( parameterNames[2], "t0");
		parametersMax[2]=1e-6;
		parameters[2]=1e-6;
		break;
	}

	/*void (*f)(int, double*, double*, int, double *) = Brix2_4Parameters_Delay;
	int parameterSize = 4;
	char parameterNames[4][512] = {"v_P",	"K_PS",		"F",		"t0"};
	double parameters[4]        = {0.001, 	0, 			1000, 		0};
	double parametersMin[4]     = {0., 		0.,     	0., 		0};
	double parametersMax[4]     = {1., 		100000., 	1000000., 	1000};*/


	/*void (*f)(int, double*, double*, int, double *) = Brix2_3Parameters;
	int parameterSize = 3;
	char parameterNames[3][512] = {"v_P","K_PS","F"};
	double parameters[3]    = {0.9, 0, 40000};
	double parametersMin[3] = {0.00001, 0, 0};
	double parametersMax[3] = {1., 0., 1000000.};*/


	/*void (*f)(int, double*, double*, int, double *) = Brix2;
	int parameterSize = 4;
	char parameterNames[4][512] = {"V_P","V_I","K_PS","F"};
	 */

	//void (*f)(int, double*, double*, int, double *) = Quadratic;
	//int parameterSize = 3;
	//void (*f)(int, double*, double*, int, double *) = Cubic;
	//void (*f)(int, double*, double*, int, double *) = Linear;
	//int parameterSize = 4;
	/*double parameters[parameterSize];
	double parametersMin[parameterSize];
	double parametersMax[parameterSize];
	for( int i=0; i<parameterSize; i++){
		parametersMin[i] = 0.00001;
		parametersMax[i] = 100000;
		parameters[i] = 0.5;
	}*/

	FILE *fp;

	chdir(dirname);

	if( pixelX>=0 && pixelY>=0){
	// FIT: SINGLE CURVE
		//stepSize = 1;
		//double gamma = fitMethodParameter;//0.0001;
		//GradientDecent(sampleSize, sampleX, sampleY, parameterSize, parameters,
		//		parametersMin, parametersMax, stepSize, gamma, Smin, 10000, f);


		double lamda = fitMethodParameter;//1;
		LevenbergMarquardt(sampleSize, sampleX, sampleY, parameterSize, parameters,
		//LevenbergMarquardt(sampleSize, samplesX[pixelX][pixelY], samplesY[pixelX][pixelY], parameterSize, parameters,
				parametersMin, parametersMax, stepSize, lamda, Smin, Imax, f);

		// OUTPUT
		double fitY[sampleSize];
		f(sampleSize, sampleX, fitY, parameterSize, parameters);
		fp = fopen("output.dat", "w+");
		for (int i = 0; i < sampleSize; i++) {
			fprintf(fp, "%e %e %e %e\n", sampleX[i], sampleY[i], fitY[i],
					//Brix2AIF( sampleX[i])
					ArterialInputFunction( sampleX[i]/60.)
			);
		}
		fclose(fp);

		for( int p=0; p<parameterSize; p++){
			fprintf(stderr, "%s = %e\n", parameterNames[p], parameters[p]);
		}
	}
	else
	{

		// FIT: PARAMETER MAP

		// statistics
		double max = 0.;
		int maxX=-1, maxY=-1;

		//double parameterMap[imageWidth][imageHeight][parameterSize];
		double ***parameterMap = (double ***) malloc(sizeof(double**) * imageWidth);
#pragma omp parallel for
		for (int x = 0; x < imageWidth; x++) {
			parameterMap[x] = (double **) malloc(sizeof(double*) * imageHeight);
			for (int y = 0; y < imageHeight; y++) {
				parameterMap[x][y] = (double *) malloc(sizeof(double*) * parameterSize);

				// init
				for (int i = 0; i < parameterSize; i++)
					parameterMap[x][y][i] = parameters[i];

				// better first guess
				/*if( inverseMethode == BRIX2DELAY){
					int i=0;
					for ( ; i < sampleSize && samplesY[x][y][i] == 0; i++) ;
					parameterMap[x][y][3] = (double)i;
				}*/
				if( inverseMethode == BRIX2DELAY || inverseMethode == BRIX2NOPERM || inverseMethode == TOFTS){
					int max_index = 0;

					// global max
					/*double max_value = 0;
					for( int i=0; i<sampleSize; i++){
						if( max_value < samplesY[x][y][i]){
							max_value = samplesY[x][y][i];
							max_index = i;
						}
					}*/

					// first non-zero value
					for( max_index=0; max_index<sampleSize && samplesY[x][y][max_index]==0; max_index++) ;

					switch( inverseMethode){
					case BRIX2DELAY:
						parameterMap[x][y][3] = (double)max_index; break;
					case BRIX2NOPERM:
						parameterMap[x][y][1] = (double)max_index; break;
					case TOFTS:
						parameterMap[x][y][2] = (double)max_index; break;
					}
				}

				//stepSize = 0.001;
				/*double gamma = fitMethodParameter;
				GradientDecent(sampleSize, samplesX[x][y], samplesY[x][y],
						parameterSize, parameterMap[x][y], parametersMin,
						parametersMax, stepSize, gamma, Smin, Imax, f);*/

				double lamda = fitMethodParameter;//1;
				LevenbergMarquardt(sampleSize, samplesX[x][y], samplesY[x][y],
						parameterSize, parameterMap[x][y], parametersMin,
						parametersMax, stepSize, lamda, Smin, Imax, f);
			}
		}

		// OUTPUT TO SEPERATE FILES
		for (int i = 0; i < parameterSize; i++) {

			char filename[512];
			sprintf(filename, "%s.dat", parameterNames[i]);
			fp = fopen(filename, "w+");

			for (int x = 0; x < imageWidth; x++) {
				for (int y = 0; y < imageHeight; y++) {
					fprintf(fp, "%i %i %e\n", x, y, parameterMap[x][y][i]);
				}
			}

			fclose( fp);

		}

		// OUTPUT TO SINGLE FILE
		{

			char filename[512];
			sprintf(filename, "%s.dat", "parameters");
			fp = fopen(filename, "w+");

			fprintf(fp, "#x,y");
			for (int i = 0; i < parameterSize; i++)
				fprintf(fp, ",%s", parameterNames[i]);
			fprintf(fp, "\n");

			for (int x = 0; x < imageWidth; x++) {
				for (int y = 0; y < imageHeight; y++) {
					if( y==0 && imageHeight>1)
						fprintf(fp, "\n");
					fprintf(fp, "%i %i", x, y);
					for (int i = 0; i < parameterSize; i++)
						fprintf(fp, " %e", parameterMap[x][y][i]);
					fprintf(fp, "\n");
				}
			}

			fclose( fp);

		}

		// OUTPUT FIT CURVES
		for (int x = 0; x < imageWidth; x++)
		for (int y = 0; y < imageHeight; y++){
			// OUTPUT
			double fitY[sampleSize];
			f(sampleSize, samplesX[x][y], fitY, parameterSize, parameterMap[x][y]);

			char filename[512];
			sprintf(filename, "fit%ix%i.dat", x,y);

			fp = fopen(filename, "w+");
			for (int i = 0; i < sampleSize; i++) {
				if( max<samplesY[x][y][i]){
					max=samplesY[x][y][i];
					maxX = x;
					maxY = y;
				}

				fprintf(fp, "%e %e %e %e\n", samplesX[x][y][i], samplesY[x][y][i], fitY[i],
						//Brix2AIF( sampleX[i])
						ArterialInputFunction( samplesX[x][y][i]/60.)
				);

			}
			fclose(fp);
		}

		fprintf ( stderr, "Maximal pixel intesity over time (%lf) was found in pixel (%i,%i)\n", max, maxX, maxY);
	}

	// FREE MEMORY
	for( int x=0; x<imageWidth; x++){
		for( int y=0; y<imageHeight; y++){
			free( samplesX[x][y]);
			free( samplesY[x][y]);
		}
		free( samplesX[x]);
		free( samplesY[x]);
	}
	free( samplesX);
	free( samplesY);

}

