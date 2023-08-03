#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP mehdi_wt_cor(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"mehdi_wt_cor", (DL_FUNC) &mehdi_wt_cor, 4},
    {NULL, NULL, 0}
};

void R_init_ilincsR(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}



// #include <R.h>
// #include <Rinternals.h>
// #include <stdlib.h> // for NULL
// #include <R_ext/Rdynload.h>
// 
// void R_init_ilincsR(DllInfo *info) {
//   R_RegisterCCallable("ilincsR", "mehdi_wtcor",  (DL_FUNC) &mehdi_wt_cor);
// }


// #include "add.h"
// #include <R_ext/Rdynload.h>
// 
// void R_init_mypackage(DllInfo *info) {
//   R_RegisterCCallable("mypackage", "add",  (DL_FUNC) &add_);
// }

