#include <stdlib.h>
#include <stdio.h>
#include <complex>
#include <cmath>
#include <cstring>
#include <climits>

#include "constants.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

/*-----------------------------------------------------------------*/

struct quad_func_params
{
    int ind;
    double a;
    double omE;
    double nu;
    complex<double> c;
    complex<double> d;
};

/*-----------------------------------------------------------------*/

extern "C"
{
void calc_velocity_integral_ (int ind, double vT, double ks, double kp, double omE, double nu, double *Es, double *Ep, double *res);
}

/*-----------------------------------------------------------------*/

double func (double x, void *params)
{
quad_func_params *P = (quad_func_params *)params;

complex<double> field_fac = (P->c + (P->d)*x)*conj((P->c + (P->d)*x));

double ind_fac;

if (P->ind == 1)
{
    ind_fac = 1.0;
}
else if (P->ind == 2)
{
    ind_fac = 1.0 + x*x;
}
else if (P->ind == 3)
{
    ind_fac = 1.0 + (1.0 + x*x)*(1.0 + x*x);
}
else
{
    ind_fac = 0.0;
    fprintf (stderr, "\nunknown index!");
}

return exp(-x*x)/(pow((P->a)*x + P->omE, 2)+(P->nu)*(P->nu))*real(field_fac)*ind_fac;
}

/*-----------------------------------------------------------------*/

void calc_velocity_integral_ (int ind, double vT, double ks, double kp, double omE, double nu, double *Es, double *Ep, double *res)
{
gsl_integration_workspace *W = gsl_integration_workspace_alloc (10000);

struct quad_func_params P;

gsl_function F;
F.function = &func;
F.params = &P;

//set structure:
P.ind = ind;
P.a = sqrt(2.0)*vT*kp;
P.omE = omE;
P.nu = nu;
P.c = (omE/ks)*(Es[0] + Es[1]*I);
P.d = sqrt(2.0)*vT*(Ep[0] + Ep[1]*I);

//start for test:
/*int i;
double x;
FILE *out;

if (!(out = fopen ("func.dat", "w")))
{
    fprintf (stderr, "\nFailed to open file %s\a\n", "func.dat");
}

for (i=0; i<10000; i++)
{
    x = -1.0+0.0002*i;
    fprintf (out, "\n%le\t%le", x, func (x, (void *)&P));
}*/
//end for test:

double err;

gsl_integration_qagi (&F, 1.0e-6, 1.0e-6, 10000, W, res, &err);
//gsl_integration_qag (&F, -10.0, 10.0, 1.0e-6, 1.0e-6, 10000, 3, W, res, &err);

gsl_integration_workspace_free (W);
}

/*-----------------------------------------------------------------*/
