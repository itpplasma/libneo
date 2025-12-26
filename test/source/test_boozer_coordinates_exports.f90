program test_boozer_coordinates_exports

    use boozer_coordinates_mod, only: ns_s_B, ns_tp_B, ns_B, n_theta_B, n_phi_B
    use boozer_coordinates_mod, only: hs_B, h_theta_B, h_phi_B
    use boozer_coordinates_mod, only: s_Bcovar_tp_B, s_Bmod_B, s_Bcovar_r_B
    use boozer_coordinates_mod, only: s_delt_delp_V, s_delt_delp_B
    use boozer_coordinates_mod, only: ns_max, derf1, derf2, derf3
    use boozer_coordinates_mod, only: use_B_r, use_del_tp_B
    use libneo_kinds, only: dp

    implicit none

    logical :: ok

    ok = .true.

    ns_s_B = 3
    ns_tp_B = 3
    ns_B = 2
    n_theta_B = 4
    n_phi_B = 5

    hs_B = 1.0_dp
    h_theta_B = 1.0_dp
    h_phi_B = 1.0_dp

    use_B_r = .true.
    use_del_tp_B = .true.

    if (.not. allocated(s_Bcovar_tp_B)) then
        allocate (s_Bcovar_tp_B(2, ns_s_B + 1, ns_B))
    end if

    if (.not. allocated(s_Bmod_B)) then
        allocate (s_Bmod_B(ns_s_B + 1, ns_tp_B + 1, ns_tp_B + 1, ns_B, &
            n_theta_B, n_phi_B))
    end if

    if (.not. allocated(s_Bcovar_r_B)) then
        allocate (s_Bcovar_r_B(ns_s_B + 1, ns_tp_B + 1, ns_tp_B + 1, ns_B, &
            n_theta_B, n_phi_B))
    end if

    if (.not. allocated(s_delt_delp_V)) then
        allocate (s_delt_delp_V(2, ns_s_B + 1, ns_tp_B + 1, ns_tp_B + 1, ns_B, &
            n_theta_B, n_phi_B))
    end if

    if (.not. allocated(s_delt_delp_B)) then
        allocate (s_delt_delp_B(2, ns_s_B + 1, ns_tp_B + 1, ns_tp_B + 1, ns_B, &
            n_theta_B, n_phi_B))
    end if

    derf1 = 0.0_dp
    derf2 = 0.0_dp
    derf3 = 0.0_dp

    if (size(derf1) /= ns_max) ok = .false.
    if (size(derf2) /= ns_max) ok = .false.
    if (size(derf3) /= ns_max) ok = .false.

    if (.not. ok) then
        error stop "ERROR boozer_coordinates_mod exports are inconsistent."
    end if

end program test_boozer_coordinates_exports
