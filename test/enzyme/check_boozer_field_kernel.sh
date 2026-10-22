#!/usr/bin/env bash
# Enzyme cross-check of the FortSym Jacobian of the Boozer field kernel.
#
# Differentiates libneo's own primal field code (compute_field_components,
# metric_tensor_symflux, metric_tensor_vmec from spline_vmec_data.f90, copied
# verbatim) with Enzyme forward mode and compares against the generated
# boozer_field_jacobian leaves. Needs LLVM flang and the Enzyme plugin of the
# same LLVM major version; it is not part of the default CTest set.
#
# Usage: test/enzyme/check_boozer_field_kernel.sh [LLVM_BIN] [ENZYME_PLUGIN]
# FLANG (default: flang on PATH) must match the LLVM version of LLVM_BIN.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
src=$(cd "$here/../.." && pwd)/src
llvm=${1:-/opt/homebrew/opt/llvm/bin}
plugin=${2:-/opt/homebrew/lib/LLVMEnzyme-22.dylib}
flang=${FLANG:-flang}
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
cd "$work"
cp "$src/generated/boozer_field_jacobian_leaves.inc" .
{
    printf 'module kern\nuse iso_c_binding\nimplicit none\n'
    printf 'integer, parameter :: dp = c_double\ncontains\n'
    for r in compute_field_components metric_tensor_symflux metric_tensor_vmec; do
        awk "/   subroutine $r\\(/,/end subroutine $r\$/" \
            "$src/vmec_support/spline_vmec_data.f90"
    done
    cat <<'FEOF'
   subroutine field_f(x, athetas, f) bind(C, name="field_f")
      real(dp), intent(in) :: x(10)
      real(dp), value :: athetas
      real(dp), intent(out) :: f(3)
      real(dp) :: sqg, bctt, bctp, bcs, bcvt, bcvp
      call compute_field_components(x(1), x(2), x(3), x(4), x(5), x(6), x(7), &
          athetas, x(10), 0.0_dp, x(8), x(9), sqg, bctt, bctp, bcs, bcvt, bcvp)
      f(1) = sqrt(bctt*bcvt + bctp*bcvp)
      f(2) = bcvt*(1.0_dp + x(8))
      f(3) = bcvp + bcvt*x(9)
   end subroutine field_f
   include 'boozer_field_jacobian_leaves.inc'
   subroutine jac_sym(x, athetas, f, j) bind(C, name="jac_sym")
      real(dp), intent(in) :: x(10)
      real(dp), value :: athetas
      real(dp), intent(out) :: f(3), j(3, 10)
      call boozer_field_jacobian(x(1), x(2), x(3), x(4), x(5), x(6), x(7), &
         x(8), x(9), x(10), athetas, &
         f(1), j(1,1), j(1,2), j(1,3), j(1,4), j(1,5), j(1,6), j(1,7), j(1,8), &
         j(1,9), j(1,10), f(2), j(2,1), j(2,2), j(2,3), j(2,4), j(2,5), j(2,6), &
         j(2,7), j(2,8), j(2,9), j(2,10), f(3), j(3,1), j(3,2), j(3,3), j(3,4), &
         j(3,5), j(3,6), j(3,7), j(3,8), j(3,9), j(3,10))
   end subroutine jac_sym
end module kern
FEOF
} > kern.f90
"$flang" -O2 -c -emit-llvm kern.f90 -o kern.bc
"$llvm/clang" -O2 -c -emit-llvm "$here/enzyme_kernel_driver.c" -o drv.bc
"$llvm/llvm-link" kern.bc drv.bc -o all.bc
"$llvm/opt" all.bc -load-pass-plugin="$plugin" -passes=enzyme -o ad.bc
"$llvm/opt" -O2 ad.bc -o opt.bc
"$flang" opt.bc -o check.x
./check.x
