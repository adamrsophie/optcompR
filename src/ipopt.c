/* Optional binding to the IPOPT C interface. Compiled against IPOPT when
 * configure finds it with pkg-config (HAVE_IPOPT); otherwise every entry
 * point reports that IPOPT is unavailable. */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <string.h>

#ifdef HAVE_IPOPT
#include <IpStdCInterface.h>

/* R closures and scratch data passed to the callbacks. */
typedef struct {
  SEXP f, grad, g, jac, hess, env;
  int error;
} ipopt_data;

/* Call fun(x, ...) in R; NULL on error. */
static SEXP call_r(ipopt_data *d, SEXP fun, SEXP args) {
  int err = 0;
  SEXP call = PROTECT(LCONS(fun, args));
  SEXP val = R_tryEval(call, d->env, &err);
  UNPROTECT(1);
  if (err) {
    d->error = 1;
    return NULL;
  }
  return val;
}

static SEXP xvec(Index n, Number *x) {
  SEXP v = PROTECT(allocVector(REALSXP, n));
  memcpy(REAL(v), x, sizeof(double) * n);
  UNPROTECT(1);
  return v;
}

/* Copy an R numeric result of length k into out. */
static Bool take(ipopt_data *d, SEXP val, Number *out, Index k) {
  if (val == NULL) return 0;
  PROTECT(val);
  if (TYPEOF(val) != REALSXP || XLENGTH(val) != k) {
    UNPROTECT(1);
    d->error = 1;
    return 0;
  }
  memcpy(out, REAL(val), sizeof(double) * k);
  UNPROTECT(1);
  return 1;
}

static Bool cb_f(Index n, Number *x, Bool new_x, Number *obj, UserDataPtr ud) {
  ipopt_data *d = (ipopt_data *) ud;
  SEXP a = PROTECT(list1(xvec(n, x)));
  Bool ok = take(d, call_r(d, d->f, a), obj, 1);
  UNPROTECT(1);
  return ok;
}

static Bool cb_grad(Index n, Number *x, Bool new_x, Number *grad,
                    UserDataPtr ud) {
  ipopt_data *d = (ipopt_data *) ud;
  SEXP a = PROTECT(list1(xvec(n, x)));
  Bool ok = take(d, call_r(d, d->grad, a), grad, n);
  UNPROTECT(1);
  return ok;
}

static Bool cb_g(Index n, Number *x, Bool new_x, Index m, Number *g,
                 UserDataPtr ud) {
  ipopt_data *d = (ipopt_data *) ud;
  if (m == 0) return 1;
  SEXP a = PROTECT(list1(xvec(n, x)));
  Bool ok = take(d, call_r(d, d->g, a), g, m);
  UNPROTECT(1);
  return ok;
}

/* Structures are dense: rows by columns for the Jacobian (column major)
 * and the lower triangle by columns for the Hessian. */
static Bool cb_jac(Index n, Number *x, Bool new_x, Index m, Index nele,
                   Index *iRow, Index *jCol, Number *values, UserDataPtr ud) {
  ipopt_data *d = (ipopt_data *) ud;
  if (values == NULL) {
    Index k = 0;
    for (Index j = 0; j < n; j++)
      for (Index i = 0; i < m; i++) {
        iRow[k] = i;
        jCol[k] = j;
        k++;
      }
    return 1;
  }
  SEXP a = PROTECT(list1(xvec(n, x)));
  Bool ok = take(d, call_r(d, d->jac, a), values, nele);
  UNPROTECT(1);
  return ok;
}

static Bool cb_hess(Index n, Number *x, Bool new_x, Number obj_factor,
                    Index m, Number *lambda, Bool new_lambda, Index nele,
                    Index *iRow, Index *jCol, Number *values, UserDataPtr ud) {
  ipopt_data *d = (ipopt_data *) ud;
  if (isNull(d->hess)) return 0;
  if (values == NULL) {
    Index k = 0;
    for (Index j = 0; j < n; j++)
      for (Index i = j; i < n; i++) {
        iRow[k] = i;
        jCol[k] = j;
        k++;
      }
    return 1;
  }
  SEXP lam = PROTECT(allocVector(REALSXP, m));
  if (m) memcpy(REAL(lam), lambda, sizeof(double) * m);
  SEXP a = PROTECT(list3(xvec(n, x), ScalarReal(obj_factor), lam));
  Bool ok = take(d, call_r(d, d->hess, a), values, nele);
  UNPROTECT(2);
  return ok;
}

SEXP optcompR_ipopt_solve(SEXP x0, SEXP xl, SEXP xu, SEXP gl, SEXP gu,
                          SEXP f, SEXP grad, SEXP g, SEXP jac, SEXP hess,
                          SEXP opts, SEXP env) {
  Index n = (Index) XLENGTH(x0), m = (Index) XLENGTH(gl);
  ipopt_data d = {f, grad, g, jac, hess, env, 0};
  int exact = !isNull(hess);
  IpoptProblem p = CreateIpoptProblem(
      n, REAL(xl), REAL(xu), m, REAL(gl), REAL(gu), n * m,
      exact ? n * (n + 1) / 2 : 0, 0, cb_f, cb_g, cb_grad, cb_jac,
      cb_hess); /* required even with a limited-memory Hessian */
  if (p == NULL) error("IPOPT could not create the problem.");
  AddIpoptNumOption(p, "nlp_lower_bound_inf", -1e300);
  AddIpoptNumOption(p, "nlp_upper_bound_inf", 1e300);
  if (!exact) AddIpoptStrOption(p, "hessian_approximation", "limited-memory");
  SEXP nm = getAttrib(opts, R_NamesSymbol);
  for (R_xlen_t k = 0; k < XLENGTH(opts); k++) {
    char *key = (char *) CHAR(STRING_ELT(nm, k));
    SEXP v = VECTOR_ELT(opts, k);
    if (isString(v)) AddIpoptStrOption(p, key, (char *) CHAR(STRING_ELT(v, 0)));
    else if (isInteger(v)) AddIpoptIntOption(p, key, INTEGER(v)[0]);
    else AddIpoptNumOption(p, key, asReal(v));
  }
  SEXP x = PROTECT(duplicate(x0));
  SEXP gv = PROTECT(allocVector(REALSXP, m));
  SEXP mg = PROTECT(allocVector(REALSXP, m));
  SEXP ml = PROTECT(allocVector(REALSXP, n));
  SEXP mu = PROTECT(allocVector(REALSXP, n));
  double obj = NA_REAL;
  int status = IpoptSolve(p, REAL(x), m ? REAL(gv) : NULL, &obj,
                          m ? REAL(mg) : NULL, REAL(ml), REAL(mu), &d);
  FreeIpoptProblem(p);
  const char *names[] = {"status", "x", "objective", "mult_g", "mult_x_L",
                         "mult_x_U", "callback_error", ""};
  SEXP out = PROTECT(mkNamed(VECSXP, names));
  SET_VECTOR_ELT(out, 0, ScalarInteger(status));
  SET_VECTOR_ELT(out, 1, x);
  SET_VECTOR_ELT(out, 2, ScalarReal(obj));
  SET_VECTOR_ELT(out, 3, mg);
  SET_VECTOR_ELT(out, 4, ml);
  SET_VECTOR_ELT(out, 5, mu);
  SET_VECTOR_ELT(out, 6, ScalarLogical(d.error));
  UNPROTECT(6);
  return out;
}

SEXP optcompR_ipopt_available(void) { return ScalarLogical(1); }

#else

SEXP optcompR_ipopt_solve(SEXP x0, SEXP xl, SEXP xu, SEXP gl, SEXP gu,
                          SEXP f, SEXP grad, SEXP g, SEXP jac, SEXP hess,
                          SEXP opts, SEXP env) {
  error("optcompR was installed without IPOPT.");
  return R_NilValue;
}

SEXP optcompR_ipopt_available(void) { return ScalarLogical(0); }

#endif

static const R_CallMethodDef call_methods[] = {
  {"optcompR_ipopt_solve", (DL_FUNC) &optcompR_ipopt_solve, 12},
  {"optcompR_ipopt_available", (DL_FUNC) &optcompR_ipopt_available, 0},
  {NULL, NULL, 0}
};

void R_init_optcompR(DllInfo *dll) {
  R_registerRoutines(dll, NULL, call_methods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
