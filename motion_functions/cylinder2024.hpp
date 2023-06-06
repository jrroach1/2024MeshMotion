// COMPILE:
//   compile:  g++ -O3 -I. -c cylinder2024.cpp
//   produces: cylinder2024.o
//  
//   cylinder2024.o can be linked with your application 
//   and cylinder2024.hpp can be included to provide the function interfaces.
// 
// USAGE EXAMPLE:
//   #include <cylinder2024.hpp>
// 
//   coords(ALPHA_SHORT, xref, yref, t, x, y)
//   coords_dot(ALPHA_SHORT, xref, yref, t, dxdx, dxdy, dxdt, dydx, dydy, dydt)
// 
// 
// FURTHER EXAMPLE AND TEST:
//   See the program "CPPtest.cpp" as an example of usage.
//
//   Build test: g++ -O3 -I. -o CPPtest CPPtest.cpp cylinder2024.cpp
// ----------------------------------------------------------------

#pragma once
#include <complex>

using namespace std;

enum AlphaFunction { ALPHA_SHORT, ALPHA_LONG }; // Choice of alpha(t) function

void coords(AlphaFunction alpha, double xref, double yref, double t, double &x, double &y);
void coords_dot(AlphaFunction alpha, double xref, double yref, double t,
                double &dxdx, double &dxdy, double &dxdt,
                double &dydx, double &dydy, double &dydt);
