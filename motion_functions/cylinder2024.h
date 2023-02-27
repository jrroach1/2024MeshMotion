// COMPILE:
//   compile:  gcc -c cylinder2024.c
//   produces: cylinder2024.o
//  
//   cylinder2024.o can be linked with your application 
//   and cylinder2024.h can be included to provice the function interfaces.
// 
// USAGE EXAMPLE:
//   #include <cyliner2024.h>
// 
//   coords(ALPHA_SHORT,xref,yref,t,x,y)
//   coords_dot(ALPHA_SHORT,xref,yref,t,dxdx,dxdy,dxdt,dydx,dydy,dydt)
// 
// 
// FURTHER EXAMPLE AND TEST:
//   See the program "Ctest.c" as an example of usage.
//
//   Build test: gcc cylinder2024.c Ctest.c
// 
// NOTE:
//   Since this implementation utilizes the complex-step method
//   to obtain derivatives of the prescribed motion functions,
//   the input reference coordinates and time must be declared
//   and input as complex variables with 0 imaginary component.
// 
// ----------------------------------------------------------------
#pragma once 

#ifdef __cplusplus
extern "C" {
#endif

// Determines use of either the short-time alpha time-activation 
// function or the long-time alpha time-activation function.
// These should be used as inputs to the coords and coords_dot 
// functions.
#define ALPHA_SHORT 1
#define ALPHA_LONG 2

void coords(int alpha,double complex xref,double complex yref, double complex t, double complex *x,double complex *y);

void coords_dot(int alpha,double complex xref,double complex yref,double complex t, double *dxdx,double *dxdy,double *dxdt,double *dydx,double *dydy,double *dydt);

#ifdef __cplusplus
}
#endif
