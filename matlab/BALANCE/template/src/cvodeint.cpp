#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <rhs.h>

#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <cvode/cvode_dense.h>
#include <sundials/sundials_dense.h>
#include <sundials/sundials_types.h>

#define Ith(v,i)     NV_Ith_S(v,i-1)       // Ith numbers components 1..NEQ
#define IJth(A,i,j)  DENSE_ELEM(A,i-1,j-1) // IJth numbers rows,cols 1..NEQ

/******************************************************************************/

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
//static int Jac(int N, realtype t,
//               N_Vector y, N_Vector fy, DlsMat J, void *user_data,
//               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/******************************************************************************/

extern "C"
{
int cvodeint_ (int *Neqp, double *x1, double *x2, double *y, double *eps);
}

/******************************************************************************/

extern "C"
{
void rhs_balance_(double *, const double *, const double *);
}

/******************************************************************************/

static void PrintFinalStats(void *cvode_mem);

/******************************************************************************/

int cvodeint_ (int *Neqp, double *x1, double *x2, double *y, double *eps)
{
//(int argc, char** argv) ehemals Argumente
	
	// ===== Systemspezifische Variablen ==================================
	//const char* filename = "profiles.dat";
	//int N = 170; // length of the time grid
	//N_Vector y0 = NULL; y0 = N_VNew_Serial(neq);

	// ====== Anfangsbedingungen ==========================================
	//Ith(y0,1) = 1.0;
	//Ith(y0,2) = 0.0;
	//y0 = N_VClone(*y);
	//y0 = *y;
	
	//double *y0, *y1;
	//y1 = NV_DATA_S(y);
	//N_Vector y0;
	//y0 = NV_DATA_S(*y);	
	//y0 = N_VClone(*y);
	//N_Vector N_VMake_Serial(long int vec_length, realtype *v_data);
	//y0 = N_VMake_Serial(neq, *y);

int neq = *Neqp;

N_Vector yv = N_VMake_Serial(neq, y);
if (!yv)
{
	fprintf(stderr, "\nerror: cvodeint: y vector allocation failed!..");
  	return 1;
}


realtype t0 = *x1;
realtype tfinal = *x2;

// ====== Genauigkeit =================================================
realtype reltol = *eps;
realtype abstol = *eps;

//N_Vector abstol = NULL;
//abstol = N_VNew_Serial(neq); 
//Ith(abstol,1) = *eps; Ith(abstol,2) = *eps;

realtype t;
int flag;

// INTEGRATIONSMETHODE -----------------------------------------------------
void *cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
if (!cvode_mem)
{
  fprintf(stderr, "\nerror: int_basis_vecs_: cvodecreate failed!..");
  return 1;
}


// INITIALISIERUNG ---------------------------------------------------------
flag = CVodeInit(cvode_mem, f, t0, yv);

// RELATIVE UND ABSOLUTE GENAUIGKEIT ---------------------------------------
flag = CVodeSStolerances(cvode_mem, reltol, abstol);

// AUFRUF DER INTEGRATIONSROUTINE ------------------------------------------
flag = CVDense(cvode_mem, neq);

flag = CVode(cvode_mem, tfinal, yv, &t, CV_NORMAL);

if (flag != CV_SUCCESS)
{
   fprintf(stderr, "\nerror: cvode failed!"); 
}

PrintFinalStats(cvode_mem);

CVodeFree(&cvode_mem);

return(0);
}

/******************************************************************************/

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{

// double *dbly, *dblyd;
// dbly = NV_DATA_S(y);
// dblyd = NV_DATA_S(ydot);

rhs_balance_ (&t, N_VGetArrayPointer(y), N_VGetArrayPointer(ydot));

return 0;
}

/******************************************************************************/

static int check_flag(void *flagvalue, char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
	      funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  return(0);
}

/******************************************************************************/

static void PrintFinalStats(void *cvode_mem)
{
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
  int flag;

  flag = CVodeGetNumSteps(cvode_mem, &nst);
  check_flag(&flag, "CVodeGetNumSteps", 1);
  flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_flag(&flag, "CVodeGetNumRhsEvals", 1);
  flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
  flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_flag(&flag, "CVodeGetNumErrTestFails", 1);
  flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
  flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);

  flag = CVDlsGetNumJacEvals(cvode_mem, &nje);
  check_flag(&flag, "CVDlsGetNumJacEvals", 1);
  flag = CVDlsGetNumRhsEvals(cvode_mem, &nfeLS);
  check_flag(&flag, "CVDlsGetNumRhsEvals", 1);

  flag = CVodeGetNumGEvals(cvode_mem, &nge);
  check_flag(&flag, "CVodeGetNumGEvals", 1);

  printf("\nFinal Statistics:\n");
  printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nfeLS = %-6ld nje = %ld\n",
	 nst, nfe, nsetups, nfeLS, nje);
  printf("nni = %-6ld ncfn = %-6ld netf = %-6ld nge = %ld\n \n",
	 nni, ncfn, netf, nge);
}

/******************************************************************************/
