/* Optional binding to the IPOPT C interface. Compiled against IPOPT when
 * configure finds it with pkg-config (HAVE_IPOPT); otherwise every entry
 * point reports that IPOPT is unavailable. */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <string.h>
#include <math.h>

#ifdef HAVE_IPOPT
#include <IpStdCInterface.h>
/* the bindings below pass double buffers; a single precision IPOPT build
   is not supported and falls back to the stub */
#ifdef IPOPT_SINGLE
#undef HAVE_IPOPT
#endif
#endif

#ifdef HAVE_IPOPT

/* IPOPT 3.14 deprecates Index, Number and Bool in favor of ipindex,
 * ipnumber and bool. */
#ifdef HAVE_IPOPT_314
typedef ipindex oc_index;
typedef ipnumber oc_number;
typedef bool oc_bool;
#else
typedef Index oc_index;
typedef Number oc_number;
typedef Bool oc_bool;
#endif

/* R closures and scratch data passed to the callbacks. */
typedef struct {
  SEXP f, grad, g, jac, hess, env;
  int error;
  char msg[512]; /* first R error message from a callback */
} ipopt_data;

/* Call fun(x, ...) in R; NULL on error. */
static SEXP call_r(ipopt_data *d, SEXP fun, SEXP args) {
  int err = 0;
  SEXP call = PROTECT(LCONS(fun, args));
  SEXP val = R_tryEvalSilent(call, d->env, &err);
  UNPROTECT(1);
  if (err) {
    /* IPOPT may recover by backtracking; keep the message for the caller */
    if (!d->error) {
      strncpy(d->msg, R_curErrorBuf(), sizeof(d->msg) - 1);
      d->msg[sizeof(d->msg) - 1] = '\0';
    }
    d->error = 1;
    return NULL;
  }
  return val;
}

static SEXP xvec(oc_index n, oc_number *x) {
  SEXP v = PROTECT(allocVector(REALSXP, n));
  memcpy(REAL(v), x, sizeof(double) * n);
  UNPROTECT(1);
  return v;
}

/* Copy an R numeric result of length k into out. */
static oc_bool take(ipopt_data *d, SEXP val, oc_number *out, oc_index k) {
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

static oc_bool cb_f(oc_index n, oc_number *x, oc_bool new_x, oc_number *obj, UserDataPtr ud) {
  ipopt_data *d = (ipopt_data *) ud;
  SEXP a = PROTECT(list1(xvec(n, x)));
  oc_bool ok = take(d, call_r(d, d->f, a), obj, 1);
  UNPROTECT(1);
  return ok;
}

static oc_bool cb_grad(oc_index n, oc_number *x, oc_bool new_x, oc_number *grad,
                    UserDataPtr ud) {
  ipopt_data *d = (ipopt_data *) ud;
  SEXP a = PROTECT(list1(xvec(n, x)));
  oc_bool ok = take(d, call_r(d, d->grad, a), grad, n);
  UNPROTECT(1);
  return ok;
}

static oc_bool cb_g(oc_index n, oc_number *x, oc_bool new_x, oc_index m, oc_number *g,
                 UserDataPtr ud) {
  ipopt_data *d = (ipopt_data *) ud;
  if (m == 0) return 1;
  SEXP a = PROTECT(list1(xvec(n, x)));
  oc_bool ok = take(d, call_r(d, d->g, a), g, m);
  UNPROTECT(1);
  return ok;
}

/* Structures are dense: rows by columns for the Jacobian (column major)
 * and the lower triangle by columns for the Hessian. */
static oc_bool cb_jac(oc_index n, oc_number *x, oc_bool new_x, oc_index m, oc_index nele,
                   oc_index *iRow, oc_index *jCol, oc_number *values, UserDataPtr ud) {
  ipopt_data *d = (ipopt_data *) ud;
  if (values == NULL) {
    oc_index k = 0;
    for (oc_index j = 0; j < n; j++)
      for (oc_index i = 0; i < m; i++) {
        iRow[k] = i;
        jCol[k] = j;
        k++;
      }
    return 1;
  }
  SEXP a = PROTECT(list1(xvec(n, x)));
  oc_bool ok = take(d, call_r(d, d->jac, a), values, nele);
  UNPROTECT(1);
  return ok;
}

static oc_bool cb_hess(oc_index n, oc_number *x, oc_bool new_x, oc_number obj_factor,
                    oc_index m, oc_number *lambda, oc_bool new_lambda, oc_index nele,
                    oc_index *iRow, oc_index *jCol, oc_number *values, UserDataPtr ud) {
  ipopt_data *d = (ipopt_data *) ud;
  if (isNull(d->hess)) return 0;
  if (values == NULL) {
    oc_index k = 0;
    for (oc_index j = 0; j < n; j++)
      for (oc_index i = j; i < n; i++) {
        iRow[k] = i;
        jCol[k] = j;
        k++;
      }
    return 1;
  }
  SEXP lam = PROTECT(allocVector(REALSXP, m));
  if (m) memcpy(REAL(lam), lambda, sizeof(double) * m);
  SEXP xv = PROTECT(xvec(n, x));
  SEXP of = PROTECT(ScalarReal(obj_factor));
  SEXP a = PROTECT(list3(xv, of, lam));
  oc_bool ok = take(d, call_r(d, d->hess, a), values, nele);
  UNPROTECT(4);
  return ok;
}

SEXP optcompR_ipopt_solve(SEXP x0, SEXP xl, SEXP xu, SEXP gl, SEXP gu,
                          SEXP f, SEXP grad, SEXP g, SEXP jac, SEXP hess,
                          SEXP opts, SEXP env) {
  oc_index n = (oc_index) XLENGTH(x0), m = (oc_index) XLENGTH(gl);
  ipopt_data d = {f, grad, g, jac, hess, env, 0, ""};
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
    oc_bool ok = 0;
    if (XLENGTH(v) != 1) {
      ok = 0;
    } else if (isString(v)) {
      ok = AddIpoptStrOption(p, key, (char *) CHAR(STRING_ELT(v, 0)));
    } else if (isInteger(v) || isReal(v)) {
      double val = asReal(v);
      /* IPOPT options are typed; try the type of the R value first, then
         the other one, since R numbers like 10 may name integer options */
      int whole = R_FINITE(val) && fabs(val) < 2e9 && val == (double) (int) val;
      if (isInteger(v)) {
        ok = AddIpoptIntOption(p, key, (oc_index) val);
        if (!ok) ok = AddIpoptNumOption(p, key, val);
      } else {
        ok = AddIpoptNumOption(p, key, val);
        if (!ok && whole) ok = AddIpoptIntOption(p, key, (oc_index) val);
      }
    }
    if (!ok) {
      FreeIpoptProblem(p);
      error("IPOPT rejected option `%s`.", key);
    }
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
                         "mult_x_U", "callback_error", "callback_message",
                         ""};
  SEXP out = PROTECT(mkNamed(VECSXP, names));
  SET_VECTOR_ELT(out, 0, ScalarInteger(status));
  SET_VECTOR_ELT(out, 1, x);
  SET_VECTOR_ELT(out, 2, ScalarReal(obj));
  SET_VECTOR_ELT(out, 3, mg);
  SET_VECTOR_ELT(out, 4, ml);
  SET_VECTOR_ELT(out, 5, mu);
  SET_VECTOR_ELT(out, 6, ScalarLogical(d.error));
  SET_VECTOR_ELT(out, 7, mkString(d.msg));
  UNPROTECT(6);
  return out;
}

/* 2: IPOPT 3.14 or later (max_wall_time exists), 1: older IPOPT */
SEXP optcompR_ipopt_available(void) {
#ifdef HAVE_IPOPT_314
  return ScalarInteger(2);
#else
  return ScalarInteger(1);
#endif
}

#else

SEXP optcompR_ipopt_solve(SEXP x0, SEXP xl, SEXP xu, SEXP gl, SEXP gu,
                          SEXP f, SEXP grad, SEXP g, SEXP jac, SEXP hess,
                          SEXP opts, SEXP env) {
  error("optcompR was installed without IPOPT.");
  return R_NilValue;
}

SEXP optcompR_ipopt_available(void) { return ScalarInteger(0); }

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
