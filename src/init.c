#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP mcmcse_inseq(SEXP, SEXP);
extern SEXP mcmcse_mbmc(SEXP, SEXP);
extern SEXP mcmcse_msvec(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"mcmcse_inseq", (DL_FUNC) &mcmcse_inseq, 2},
    {"mcmcse_mbmc",  (DL_FUNC) &mcmcse_mbmc,  2},
    {"mcmcse_msvec", (DL_FUNC) &mcmcse_msvec, 3},
    {NULL, NULL, 0}
};

void R_init_mcmcse(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
