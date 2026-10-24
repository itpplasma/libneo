/* Enzyme forward mode of libneo's field kernel versus the FortSym Jacobian. */
#include <stdio.h>
#include <math.h>
extern void field_f(const double *x, double at, double *f);
extern void jac_sym(const double *x, double at, double *f, double *j);
extern int enzyme_dup, enzyme_const, enzyme_dupnoneed;
void __enzyme_fwddiff(void *, ...);
int main(void) {
  /* point with libneo magnitudes (cm, G cm^2) */
  double x[10] = {1000.0, 30.0, 5.0, -40.0, 25.0, 95.0, 12.0, 0.08, -0.03, -2.1e9};
  double at = -5.0e9, f[3], fs[3], js[30], maxrel = 0;
  jac_sym(x, at, fs, js);
  for (int k = 0; k < 10; ++k) {
    double dx[10] = {0}, df[3] = {0};
    dx[k] = 1.0;
    __enzyme_fwddiff((void *)field_f, enzyme_dup, x, dx, enzyme_const, at, enzyme_dup, f, df);
    for (int i = 0; i < 3; ++i) {
      double ref = js[i + 3 * k], scale = 0;
      for (int kk = 0; kk < 10; ++kk) scale = fmax(scale, fabs(js[i + 3 * kk]));
      double r = fabs(df[i] - ref) / scale;
      if (r > maxrel) maxrel = r;
    }
  }
  printf("f(enzyme primal) %.17g %.17g %.17g\n", f[0], f[1], f[2]);
  printf("f(fortsym)       %.17g %.17g %.17g\n", fs[0], fs[1], fs[2]);
  printf("max row-scaled |J_enzyme - J_fortsym| = %.3e\n", maxrel);
  return maxrel < 1e-13 ? 0 : 1;
}
