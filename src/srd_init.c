#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(srd)(int *r, int *cd, double *loci);

static const R_FortranMethodDef FortranEntries[] = {
    {"srd", (DL_FUNC) &F77_NAME(srd), 14},
    {NULL, NULL, 0}
};

void R_init_srd(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
