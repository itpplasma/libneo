#!/usr/bin/env python3
"""
Systematic verification of Cerfon-Freidberg basis function derivatives

Compares analytical_gs_solver.f90 implementation against symbolic computation.
This script checks EVERY SINGLE derivative term-by-term.

References:
- Cerfon & Freidberg, Physics of Plasmas 17, 032502 (2010)
- Verena Eslbauer, Bachelor thesis, TU Graz (2017)
"""

import sympy as sp
from sympy import symbols, diff, log, simplify, expand

# Symbolic variables
x, y = symbols('x y', real=True, positive=True)
A = symbols('A', real=True)

# Define the 7 basis functions (Cerfon-Freidberg Eq. 8)
BASIS_FUNCTIONS = {
    1: ('1', 1),
    2: ('x²', x**2),
    3: ('y² - x²·ln(x)', y**2 - x**2 * log(x)),
    4: ('x⁴ - 4x²y²', x**4 - 4*x**2*y**2),
    5: ('2y⁴ - 9x²y² + [3x⁴ - 12x²y²]·ln(x)',
        2*y**4 - 9*x**2*y**2 + (3*x**4 - 12*x**2*y**2)*log(x)),
    6: ('x⁶ - 12x⁴y² + 8x²y⁴',
        x**6 - 12*x**4*y**2 + 8*x**2*y**4),
    7: ('8y⁶ - 140x²y⁴ + 75x⁴y² + [180x⁴y² - 120x²y⁴ - 15x⁶]·ln(x)',
        8*y**6 - 140*x**2*y**4 + 75*x**4*y**2 + (180*x**4*y**2 - 120*x**2*y**4 - 15*x**6)*log(x))
}

# Particular solution (Cerfon-Freidberg Eq. 9)
PSI_P = ('x⁴/8 + A·[x²·ln(x)/2 - x⁴/8]', x**4/8 + A*(x**2*log(x)/2 - x**4/8))

# Fortran implementation (from analytical_gs_solver.f90)
FORTRAN_IMPL = {
    # First derivatives dpsi_i/dx
    'dpsi_1_dx': '0.0_dp',
    'dpsi_2_dx': '2.0_dp * x',
    'dpsi_3_dx': '-2.0_dp * x * log(x) - x',
    'dpsi_4_dx': '4.0_dp * x**3 - 8.0_dp * x * y**2',
    'dpsi_5_dx': ('-18.0_dp * x * y**2 + 12.0_dp * x**3 * logx + 3.0_dp * x**3 '
                  '- 24.0_dp * x * y**2 * logx - 12.0_dp * x * y**2'),
    'dpsi_6_dx': '6.0_dp * x**5 - 48.0_dp * x**3 * y**2 + 16.0_dp * x * y**4',
    'dpsi_7_dx': ('-280.0_dp * x * y**4 + 150.0_dp * x**3 * y**2 '
                  '+ 720.0_dp * x**3 * y**2 * logx + 180.0_dp * x**3 * y**2 '
                  '- 240.0_dp * x * y**4 * logx - 120.0_dp * x * y**4 '
                  '- 90.0_dp * x**5 * logx - 15.0_dp * x**5'),

    # First derivatives dpsi_i/dy
    'dpsi_1_dy': '0.0_dp',
    'dpsi_2_dy': '0.0_dp',
    'dpsi_3_dy': '2.0_dp * y',
    'dpsi_4_dy': '-8.0_dp * x**2 * y',
    'dpsi_5_dy': '8.0_dp * y**3 - 18.0_dp * x**2 * y - 24.0_dp * x**2 * y * logx',
    'dpsi_6_dy': '-24.0_dp * x**4 * y + 32.0_dp * x**2 * y**3',
    'dpsi_7_dy': ('48.0_dp * y**5 - 560.0_dp * x**2 * y**3 + 150.0_dp * x**4 * y '
                  '+ 360.0_dp * x**4 * y * logx - 480.0_dp * x**2 * y**3 * logx'),

    # Second derivatives d2psi_i/dx²
    'd2psi_1_dxx': '0.0_dp',
    'd2psi_2_dxx': '2.0_dp',
    'd2psi_3_dxx': '-2.0_dp * log(x) - 3.0_dp',
    'd2psi_4_dxx': '12.0_dp * x**2 - 8.0_dp * y**2',
    'd2psi_5_dxx': '-18.0_dp * y**2 + 36.0_dp * x**2 * logx + 21.0_dp * x**2 - 24.0_dp * y**2 * logx - 36.0_dp * y**2',
    'd2psi_6_dxx': '30.0_dp * x**4 - 144.0_dp * x**2 * y**2 + 16.0_dp * y**4',
    'd2psi_7_dxx': ('-280.0_dp * y**4 + 450.0_dp * x**2 * y**2 '
                    '+ 2160.0_dp * x**2 * y**2 * logx + 1080.0_dp * x**2 * y**2 '
                    '- 240.0_dp * y**4 * logx - 360.0_dp * y**4 '
                    '- 450.0_dp * x**4 * logx - 165.0_dp * x**4'),

    # Second derivatives d2psi_i/dy²
    'd2psi_1_dyy': '0.0_dp',
    'd2psi_2_dyy': '0.0_dp',
    'd2psi_3_dyy': '2.0_dp',
    'd2psi_4_dyy': '-8.0_dp * x**2',
    'd2psi_5_dyy': '24.0_dp * y**2 - 18.0_dp * x**2 - 24.0_dp * x**2 * logx',
    'd2psi_6_dyy': '-24.0_dp * x**4 + 96.0_dp * x**2 * y**2',
    'd2psi_7_dyy': ('240.0_dp * y**4 - 1680.0_dp * x**2 * y**2 + 150.0_dp * x**4 '
                    '+ 360.0_dp * x**4 * logx - 1440.0_dp * x**2 * y**2 * logx'),

    # Particular solution derivatives
    'dpsi_p_dx': 'x**3 / 2.0_dp + A_param * (x * log(x) + x / 2.0_dp - x**3 / 2.0_dp)',
    'dpsi_p_dy': '0.0_dp',
    'd2psi_p_dxx': '1.5_dp * x**2 + A_param * (log(x) + 1.5_dp - 1.5_dp * x**2)',
    'd2psi_p_dyy': '0.0_dp',
}

def fortran_to_sympy(fortran_str):
    """Convert Fortran expression string to SymPy expression"""
    # Replace Fortran syntax with Python/SymPy
    expr_str = fortran_str.replace('_dp', '')
    expr_str = expr_str.replace('logx', 'log(x)')
    expr_str = expr_str.replace('A_param', 'A')
    try:
        return sp.sympify(expr_str)
    except:
        print(f"ERROR parsing: {fortran_str}")
        return None

def check_derivative(name, symbolic, fortran_str, test_points=None):
    """Check if Fortran implementation matches symbolic derivative"""
    print(f"\nChecking {name}")
    print("="*70)

    fortran_expr = fortran_to_sympy(fortran_str)
    if fortran_expr is None:
        print("❌ FAILED TO PARSE FORTRAN")
        return False

    # Remove benign floating-point noise (2.0 -> 2) before comparison
    fortran_expr = sp.nsimplify(fortran_expr)

    symbolic_expanded = expand(symbolic)
    fortran_expanded = expand(fortran_expr)
    symbolic_expanded = sp.nsimplify(symbolic_expanded)
    fortran_expanded = sp.nsimplify(fortran_expanded)

    print(f"Symbolic: {symbolic_expanded}")
    print(f"Fortran:  {fortran_expanded}")

    diff_expr = simplify(symbolic_expanded - fortran_expanded)
    diff_expr = sp.nsimplify(diff_expr)

    if diff_expr == 0 or diff_expr.equals(0):
        match = True
        diff_expr = sp.Integer(0)
    else:
        diff_expr_simplified = sp.simplify(sp.expand(diff_expr))
        if diff_expr_simplified == 0 or diff_expr_simplified.equals(0):
            match = True
            diff_expr = sp.Integer(0)
        else:
            subs_map = {x: sp.Rational(3, 2), y: sp.Rational(2, 5)}
            if A in diff_expr_simplified.free_symbols:
                subs_map[A] = sp.Rational(-1, 7)
            try:
                test_val = sp.simplify(diff_expr_simplified.subs(subs_map))
            except Exception:  # noqa: BLE001
                test_val = None

            if test_val == 0:
                match = True
                diff_expr = sp.Integer(0)
            else:
                try:
                    numeric_expr = sp.N(diff_expr_simplified.subs(subs_map), 50)
                except Exception:  # noqa: BLE001
                    numeric_expr = None

                if numeric_expr is not None:
                    if getattr(numeric_expr, "is_zero", False):
                        match = True
                        diff_expr = sp.Integer(0)
                        numeric_val = 0.0
                    else:
                        try:
                            numeric_val = float(abs(numeric_expr.evalf()))
                        except Exception:  # noqa: BLE001
                            numeric_val = None
                else:
                    numeric_val = None

                if numeric_val is not None and numeric_val < 1.0e-40:
                    match = True
                    diff_expr = sp.Integer(0)
                else:
                    diff_expr = diff_expr_simplified
                    match = False

    if match:
        print("✅ MATCH")
    else:
        print(f"❌ MISMATCH")
        print(f"Difference: {diff_expr}")

        # Try numerical evaluation
        if test_points:
            print("\nNumerical check at test points:")
            for x_val, y_val in test_points:
                try:
                    sym_val = float(symbolic.subs([(x, x_val), (y, y_val)]))
                    fort_val = float(fortran_expr.subs([(x, x_val), (y, y_val)]))
                    print(f"  ({x_val}, {y_val}): sym={sym_val:.12f}, fort={fort_val:.12f}, diff={abs(sym_val-fort_val):.2e}")
                except:
                    print(f"  ({x_val}, {y_val}): Could not evaluate numerically")

    return match

def main():
    print("="*80)
    print("SYSTEMATIC VERIFICATION OF CERFON-FREIDBERG DERIVATIVES")
    print("="*80)
    print("\nComparing analytical_gs_solver.f90 against SymPy symbolic computation")
    print()

    test_points = [(1.2, 0.3), (1.0, 0.0), (1.5, 0.5)]
    all_passed = True
    results = {}

    # Check each basis function
    for i in range(1, 8):
        desc, psi_i = BASIS_FUNCTIONS[i]
        print(f"\n{'='*80}")
        print(f"ψ_{i} = {desc}")
        print('='*80)

        results[f'psi_{i}'] = {}

        # First derivative wrt x
        dpsi_dx_sym = diff(psi_i, x)
        passed = check_derivative(
            f'dpsi_{i}_dx',
            dpsi_dx_sym,
            FORTRAN_IMPL[f'dpsi_{i}_dx'],
            test_points
        )
        results[f'psi_{i}']['dpsi_dx'] = passed
        all_passed = all_passed and passed

        # First derivative wrt y
        dpsi_dy_sym = diff(psi_i, y)
        passed = check_derivative(
            f'dpsi_{i}_dy',
            dpsi_dy_sym,
            FORTRAN_IMPL[f'dpsi_{i}_dy'],
            test_points
        )
        results[f'psi_{i}']['dpsi_dy'] = passed
        all_passed = all_passed and passed

        # Second derivative wrt x
        d2psi_dxx_sym = diff(psi_i, x, 2)
        passed = check_derivative(
            f'd2psi_{i}_dxx',
            d2psi_dxx_sym,
            FORTRAN_IMPL[f'd2psi_{i}_dxx'],
            test_points
        )
        results[f'psi_{i}']['d2psi_dxx'] = passed
        all_passed = all_passed and passed

        # Second derivative wrt y
        d2psi_dyy_sym = diff(psi_i, y, 2)
        passed = check_derivative(
            f'd2psi_{i}_dyy',
            d2psi_dyy_sym,
            FORTRAN_IMPL[f'd2psi_{i}_dyy'],
            test_points
        )
        results[f'psi_{i}']['d2psi_dyy'] = passed
        all_passed = all_passed and passed

    # Check particular solution
    print(f"\n{'='*80}")
    print(f"ψ_p (Particular Solution)")
    print('='*80)

    desc, psi_p = PSI_P
    results['psi_p'] = {}

    # Need to include A in test points
    test_points_A = [(1.2, 0.3, -0.155), (1.0, 0.0, -0.155)]

    dpsi_p_dx_sym = diff(psi_p, x)
    passed = check_derivative(
        'dpsi_p_dx',
        dpsi_p_dx_sym,
        FORTRAN_IMPL['dpsi_p_dx']
    )
    results['psi_p']['dpsi_dx'] = passed
    all_passed = all_passed and passed

    dpsi_p_dy_sym = diff(psi_p, y)
    passed = check_derivative(
        'dpsi_p_dy',
        dpsi_p_dy_sym,
        FORTRAN_IMPL['dpsi_p_dy']
    )
    results['psi_p']['dpsi_dy'] = passed
    all_passed = all_passed and passed

    d2psi_p_dxx_sym = diff(psi_p, x, 2)
    passed = check_derivative(
        'd2psi_p_dxx',
        d2psi_p_dxx_sym,
        FORTRAN_IMPL['d2psi_p_dxx']
    )
    results['psi_p']['d2psi_dxx'] = passed
    all_passed = all_passed and passed

    d2psi_p_dyy_sym = diff(psi_p, y, 2)
    passed = check_derivative(
        'd2psi_p_dyy',
        d2psi_p_dyy_sym,
        FORTRAN_IMPL['d2psi_p_dyy']
    )
    results['psi_p']['d2psi_dyy'] = passed
    all_passed = all_passed and passed

    # Summary
    print("\n" + "="*80)
    print("VERIFICATION SUMMARY")
    print("="*80)

    for func_name, derivatives in results.items():
        print(f"\n{func_name}:")
        for deriv_name, passed in derivatives.items():
            status = "✅" if passed else "❌"
            print(f"  {status} {deriv_name}")

    print("\n" + "="*80)
    if all_passed:
        print("✅ ALL DERIVATIVES VERIFIED SUCCESSFULLY")
    else:
        print("❌ SOME DERIVATIVES HAVE ERRORS - SEE DETAILS ABOVE")
    print("="*80)

    return 0 if all_passed else 1

if __name__ == '__main__':
    import sys
    sys.exit(main())
