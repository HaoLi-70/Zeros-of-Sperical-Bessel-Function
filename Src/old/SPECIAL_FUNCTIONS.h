
#ifndef SPECIAL_FUNCTIONS_h
#define SPECIAL_FUNCTIONS_h

/******************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "ALLOCATION.h"

/******************************************************************************************/


#define M_MAX(x,y) ({      \
typeof(x) _x=(x);    \
typeof(y) _y=(y);    \
(void) (&_x==&_y);   \
_x>_y?_x:_y;         \
})

#define M_MIN(x,y) ({      \
typeof(x) _x=(x);    \
typeof(y) _y=(y);    \
(void) (&_x==&_y);   \
_x>_y?_y:_x;         \
})

#define M_SWAP(x,y) ({     \
(void) (&x==&y);     \
typeof(x) _x;        \
_x=(x);              \
x=(y);               \
y=_x;                \
})

/******************************************************************************************/

extern double Factorial(int n);

extern double Factorial_DB(int n);

extern double Legendre_Hoeksema(int L, int M, double Theta);

extern double Legendre_Hoeksema_P(int L, int M, double Theta);

extern double Associated_Legendre(int l, int m, double x);

extern double Harmonic_Coefficient(int L, int M);

extern complex double Spherical_Harmonic(int L, int M, double X, double Phi);

extern complex double Spherical_Harmonic_conjugate(int L, int M, double X, double Phi);

extern void Spherical_Bessel(double x, int Lmax, double *J, double *Jp);

extern void Compute_SB_Zeros(double **Zeros_Spherical_Bessel, int Zero_L, int Zero_N, \
                             int Flag_Derivative, int Flag_X01);

extern void Sph_Bessel_Output(char *filename, double x0, double x1, int num, int indx_l);

extern void Zeros_Output(char *filename, int Zero_L, int Zero_N, int Flag_Derivative, \
                         int Flag_X01);

/******************************************************************************************/

#endif /* SPECIAL_FUNCTIONS_h */
