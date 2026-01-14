module neo_bspline
    !! Umbrella module for B-spline functionality.
    !! Re-exports types and procedures from submodules:
    !!   - neo_bspline_base: Core types, init, eval, basis functions
    !!   - neo_bspline_interp: Direct interpolation (collocation solve)
    !!   - neo_bspline_lsq: Least-squares fitting (CGLS solver)
    use neo_bspline_base
    use neo_bspline_interp
    use neo_bspline_lsq
    implicit none
end module neo_bspline
