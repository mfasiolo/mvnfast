#include <Rinternals.h>
#include <R_ext/Rdynload.h>

SEXP dmixtCpp(SEXP X_, SEXP mu_, SEXP sigma_, SEXP df_, SEXP w_, SEXP log_, SEXP ncores_, SEXP isChol_, SEXP A_); 

SEXP dmvtCpp(SEXP X_, SEXP mu_, SEXP sigma_, SEXP df_, SEXP log_, SEXP ncores_, SEXP isChol_);

SEXP mahaCpp(SEXP X, SEXP mu, SEXP sigma, SEXP ncores, SEXP isChol);

SEXP msCpp(SEXP init_, SEXP X_, SEXP cholDec_, SEXP ncores_, SEXP tol_, SEXP store_);

SEXP rmixnCpp(SEXP n_, SEXP mu_, SEXP sigma_, SEXP indV_, SEXP ncores_, SEXP isChol_,  SEXP retInd_, SEXP A_);

SEXP rmixtCpp(SEXP n_, SEXP mu_, SEXP sigma_, SEXP df_, SEXP indV_, SEXP ncores_, SEXP isChol_,  SEXP retInd_, SEXP A_);

SEXP rmvnCpp(SEXP n_, SEXP mu_, SEXP sigma_, SEXP ncores_, SEXP isChol_, SEXP A_); 

SEXP rmvtCpp(SEXP n_, SEXP mu_, SEXP sigma_, SEXP df_, SEXP ncores_, SEXP isChol_, SEXP A_);

static const R_CallMethodDef CallEntries[] = {
  {"dmixtCpp", (DL_FUNC) &dmixtCpp, 9},
  {"dmvtCpp", (DL_FUNC) &dmvtCpp, 7},
  {"mahaCpp", (DL_FUNC) &mahaCpp, 5},
  {"msCpp", (DL_FUNC) &msCpp, 6},
  {"rmixnCpp", (DL_FUNC) &rmixnCpp, 8},
  {"rmixtCpp", (DL_FUNC) &rmixtCpp, 9},
  {"rmvnCpp", (DL_FUNC) &rmvnCpp, 6},
  {"rmvtCpp", (DL_FUNC) &rmvtCpp, 7},
  {NULL, NULL, 0}
};

void R_init_mvnfast(DllInfo *info)
{
  R_registerRoutines(info, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}



