module efit_to_boozer

   use efit_to_boozer_mod, only: psimax, nlabel, ntheta, nspl, twopi, rbeg, &
                                 rsmall, qsaf, psi_pol, psi_tor_vac, psi_tor, C_const, R_spl, Z_spl, &
                                 bmod_spl, sqgnorm_spl, Gfunc_spl, rmn, rmx, zmn, zmx, raxis, zaxis
   use field_mod, only: icall
   use field_c_mod, only: icall_c
   use field_eq_mod, only: reset_field_eq_state

   implicit none

   double precision, dimension(:), allocatable :: R_oft, Z_oft, al_oft, B_oft
   double precision, dimension(:), allocatable :: Rmn_c, Rmn_s, Zmn_c, Zmn_s, almn_c, almn_s, Bmn_c, Bmn_s
   double complex, dimension(:), allocatable :: calE, calEm_times_Jb
   logical :: initialized = .false.

contains

   subroutine init(construction, method, allow_poor_equilibrium, override_reason, exp_bpol)
      integer, intent(in), optional :: construction, allow_poor_equilibrium, exp_bpol
      !f2py integer optional, intent(in) :: construction = 0
      !f2py integer optional, intent(in) :: allow_poor_equilibrium = -1
      !f2py integer optional, intent(in) :: exp_bpol = -1
      character(len=*), intent(in), optional :: method, override_reason
      integer :: code
      character(len=4096) :: message

      call init_contour(construction, method, allow_poor_equilibrium, override_reason, exp_bpol, code, message)
      if (code /= 0) then
         call deinit
         call python_error(trim(message))
      else
         initialized = .true.
      end if
   end subroutine init

   subroutine init_contour(construction, method, allow_poor_equilibrium, override_reason, exp_bpol, code, message)
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      use efit_sfl, only: efit_sfl_t, efit_sfl_config_t, efit_sfl_status_t, efit_sfl_init, EFIT_SFL_OK
      use efit_sfl_magdata, only: efit_to_boozer_read_options, efit_sfl_magdata_fill
      integer, intent(in), optional :: construction, allow_poor_equilibrium, exp_bpol
      character(len=*), intent(in), optional :: method, override_reason
      integer, intent(out) :: code
      character(len=4096), intent(out) :: message
      type(efit_sfl_t) :: equilibrium
      type(efit_sfl_config_t) :: cfg
      type(efit_sfl_status_t) :: status
      integer :: unit, ios, nstep, nsurfmax, nsurf, mpol, j, construction_inp, write_gs_map
      character(len=1024) :: gfile, line

      call deinit
      code = 7
      message = 'Cannot read the seven-record efit_to_boozer.inp'
      open(newunit=unit, file='efit_to_boozer.inp', status='old', action='read', iostat=ios)
      if (ios /= 0) return
      read(unit, *, iostat=ios) nstep
      if (ios == 0) read(unit, *, iostat=ios) nlabel
      if (ios == 0) read(unit, *, iostat=ios) ntheta
      if (ios == 0) read(unit, *, iostat=ios) nsurfmax
      if (ios == 0) read(unit, *, iostat=ios) nsurf
      if (ios == 0) read(unit, *, iostat=ios) mpol
      if (ios == 0) read(unit, *, iostat=ios) psimax
      if (ios == 0) then
         cfg%n_theta = max(128, ntheta + modulo(ntheta, 2))
         call efit_to_boozer_read_options(unit, construction_inp, cfg, write_gs_map)
      end if
      close(unit)
      if (ios /= 0) return
      if (.not. ieee_is_finite(psimax)) return
      message = 'Invalid EFIT grid sizes in efit_to_boozer.inp'
      if (nlabel < 4 .or. ntheta < 16 .or. nsurf < 1 .or. mpol < 0 .or. &
          nstep < 1 .or. nsurfmax < 1) return
      if (present(construction)) then
         if (construction /= 0) construction_inp = construction
      end if
      message = 'EFIT construction must be 0 (deck), 1 (contour), or 2 (legacy)'
      if (construction_inp < 1 .or. construction_inp > 2) return
      if (present(method)) then
         select case (trim(method))
         case ('contour')
            construction_inp = 1
         case ('legacy')
            construction_inp = 2
         case ('')
         case default
            message = 'EFIT method must be contour or legacy'
            return
         end select
      end if
      if (present(exp_bpol)) then
         message = 'exp_bpol must be -1 (infer), 0 (Wb/rad), or 1 (Wb)'
         if (exp_bpol < -1 .or. exp_bpol > 1) return
         if (exp_bpol >= 0) cfg%exp_bpol = exp_bpol
      end if
      message = 'Cannot read the g-file path from field_divB0.inp record 7'
      open(newunit=unit, file='field_divB0.inp', status='old', action='read', iostat=ios)
      if (ios /= 0) return
      do j = 1, 6
         read(unit, '(a)', iostat=ios) line
         if (ios /= 0) exit
      end do
      if (ios == 0) read(unit, *, iostat=ios) gfile
      close(unit)
      if (ios /= 0) return
      if (construction_inp == 2) then
         call init_legacy
         code = EFIT_SFL_OK
         message = ''
         return
      end if
      if (psimax /= 1.0d10 .and. psimax < 1.0d299) cfg%psi_edge_rel = psimax
      if (present(allow_poor_equilibrium)) then
         if (allow_poor_equilibrium >= 0) cfg%allow_poor_equilibrium = allow_poor_equilibrium /= 0
      end if
      if (present(override_reason)) then
         if (len_trim(override_reason) > 0) cfg%override_reason = override_reason
      end if
      if (cfg%allow_poor_equilibrium .and. len_trim(cfg%override_reason) == 0) then
         message = 'allow_poor_equilibrium requires an explicit override_reason'
         return
      end if
      call efit_sfl_init(equilibrium, trim(gfile), cfg, status)
      code = status%code
      if (allocated(status%message)) message = status%message
      if (code /= EFIT_SFL_OK .and. .not. status%overridden) return
      allocate(rbeg(nlabel), rsmall(nlabel), qsaf(nlabel), psi_pol(0:nlabel))
      allocate(psi_tor_vac(nlabel), psi_tor(0:nlabel), C_const(nlabel))
      allocate(R_spl(0:nspl, 0:ntheta, nlabel), Z_spl(0:nspl, 0:ntheta, nlabel), &
         bmod_spl(0:nspl, 0:ntheta, nlabel), sqgnorm_spl(0:nspl, 0:ntheta, nlabel), &
         Gfunc_spl(0:nspl, 0:ntheta, nlabel))
      ! Accepted quality diagnostics stay in equilibrium%diag; grid work starts clean.
      status = efit_sfl_status_t()
      call efit_sfl_magdata_fill(equilibrium, status)
      if (status%code /= EFIT_SFL_OK) then
         code = status%code
         if (allocated(status%message)) message = status%message
         return
      end if
      call spline_magdata_in_symfluxcoord
      allocate(calE(ntheta), calEm_times_Jb(ntheta))
      allocate(R_oft(ntheta), Z_oft(ntheta), al_oft(ntheta), B_oft(ntheta))
      allocate(Rmn_c(0:mpol), Rmn_s(0:mpol), Zmn_c(0:mpol), Zmn_s(0:mpol))
      allocate(almn_c(0:mpol), almn_s(0:mpol), Bmn_c(0:mpol), Bmn_s(0:mpol))
      code = EFIT_SFL_OK
      message = ''
   end subroutine init_contour

   subroutine python_error(message)
      use, intrinsic :: iso_c_binding, only: c_char, c_null_char
      character(len=*), intent(in) :: message
      interface
         subroutine set_python_error(text) bind(C, name='libneo_efit_python_error')
            import c_char
            character(kind=c_char), intent(in) :: text(*)
         end subroutine set_python_error
      end interface
      call set_python_error(trim(message)//c_null_char)
   end subroutine python_error

   subroutine contour_transform(gfile, s_list, mpol, n_theta, n_q, route, &
      use_psi_edge, psi_edge_rel, allow_poor_equilibrium, override_reason, &
      coefficients, fluxes, state, status_code, overridden, message, exp_bpol)
      character(len=*), intent(in) :: gfile, override_reason
      double precision, intent(in) :: s_list(:), psi_edge_rel
      integer, intent(in) :: mpol, n_theta, n_q, route
      integer, intent(in), optional :: exp_bpol
      !f2py integer optional, intent(in) :: exp_bpol = -1
      logical, intent(in) :: use_psi_edge, allow_poor_equilibrium
      double precision, intent(out) :: coefficients(8, 0:mpol, size(s_list)), &
         fluxes(3, size(s_list)), state(4)
      integer, intent(out) :: status_code
      logical, intent(out) :: overridden
      character(len=4096), intent(out) :: message
      double precision :: physical(2, size(s_list)), rminor

      call contour_transform_worker(gfile, s_list, mpol, n_theta, n_q, route, &
         use_psi_edge, psi_edge_rel, allow_poor_equilibrium, override_reason, &
         coefficients, fluxes, state, physical, rminor, status_code, overridden, message, .false., exp_bpol)
   end subroutine contour_transform

   subroutine contour_transform_extended(gfile, s_list, mpol, n_theta, n_q, route, &
      use_psi_edge, psi_edge_rel, allow_poor_equilibrium, override_reason, &
      coefficients, fluxes, state, physical, rminor, status_code, overridden, message, exp_bpol)
      character(len=*), intent(in) :: gfile, override_reason
      double precision, intent(in) :: s_list(:), psi_edge_rel
      integer, intent(in) :: mpol, n_theta, n_q, route
      integer, intent(in), optional :: exp_bpol
      !f2py integer optional, intent(in) :: exp_bpol = -1
      logical, intent(in) :: use_psi_edge, allow_poor_equilibrium
      double precision, intent(out) :: coefficients(8, 0:mpol, size(s_list)), &
         fluxes(3, size(s_list)), state(4), physical(2, size(s_list)), rminor
      integer, intent(out) :: status_code
      logical, intent(out) :: overridden
      character(len=4096), intent(out) :: message

      call contour_transform_worker(gfile, s_list, mpol, n_theta, n_q, route, &
         use_psi_edge, psi_edge_rel, allow_poor_equilibrium, override_reason, &
         coefficients, fluxes, state, physical, rminor, status_code, overridden, message, .true., exp_bpol)
   end subroutine contour_transform_extended

   subroutine contour_transform_worker(gfile, s_list, mpol, n_theta, n_q, route, &
      use_psi_edge, psi_edge_rel, allow_poor_equilibrium, override_reason, &
      coefficients, fluxes, state, physical, rminor, status_code, overridden, message, check_finite, exp_bpol)
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      use efit_sfl, only: efit_sfl_t, efit_sfl_config_t, efit_sfl_status_t, &
         efit_sfl_surface_t, efit_sfl_init, efit_sfl_transform, EFIT_SFL_OK
      character(len=*), intent(in) :: gfile, override_reason
      double precision, intent(in) :: s_list(:), psi_edge_rel
      integer, intent(in) :: mpol, n_theta, n_q, route
      integer, intent(in), optional :: exp_bpol
      !f2py integer optional, intent(in) :: exp_bpol = -1
      logical, intent(in) :: use_psi_edge, allow_poor_equilibrium, check_finite
      double precision, intent(out) :: coefficients(8, 0:mpol, size(s_list)), &
         fluxes(3, size(s_list)), state(4), physical(2, size(s_list)), rminor
      integer, intent(out) :: status_code
      logical, intent(out) :: overridden
      character(len=4096), intent(out) :: message
      type(efit_sfl_t) :: equilibrium
      type(efit_sfl_config_t) :: cfg
      type(efit_sfl_status_t) :: status
      type(efit_sfl_surface_t) :: surfaces(size(s_list))
      integer :: k

      coefficients = 0.0d0
      fluxes = 0.0d0
      state = 0.0d0
      physical = 0.0d0
      rminor = 0.0d0
      status_code = 7
      overridden = .false.
      message = 'Invalid contour transform sizes, route, or flux labels'
      if (mpol < 0 .or. n_q < 6 .or. n_theta < 16 .or. modulo(n_theta, 2) /= 0 &
          .or. route < 1 .or. route > 5 .or. size(s_list) < 1) then
         call python_error(trim(message))
         return
      end if
      if (any(.not. ieee_is_finite(s_list)) .or. &
          any(s_list <= 0.0d0) .or. any(s_list > 1.0d0)) then
         call python_error(trim(message))
         return
      end if
      if (use_psi_edge) then
         if (.not. ieee_is_finite(psi_edge_rel)) then
            call python_error('EFIT psi_edge_rel must be finite')
            return
         end if
      end if
      if (allow_poor_equilibrium .and. len_trim(override_reason) == 0) then
         message = 'allow_poor_equilibrium requires an explicit override_reason'
         call python_error(trim(message))
         return
      end if
      if (present(exp_bpol)) cfg%exp_bpol = exp_bpol
      cfg%n_theta = n_theta
      cfg%n_q = n_q
      cfg%route = route
      cfg%allow_poor_equilibrium = allow_poor_equilibrium
      cfg%override_reason = override_reason
      if (use_psi_edge) cfg%psi_edge_rel = psi_edge_rel
      call efit_sfl_init(equilibrium, gfile, cfg, status)
      status_code = status%code
      overridden = status%overridden
      if (allocated(status%message)) message = status%message
      if (status_code /= EFIT_SFL_OK .and. .not. overridden) then
         call python_error(trim(message))
         return
      end if
      call efit_sfl_transform(equilibrium, s_list, mpol, surfaces, status)
      if (status%code /= EFIT_SFL_OK) then
         status_code = status%code
         overridden = .false.
         if (allocated(status%message)) message = status%message
         call python_error(trim(message))
         return
      end if
      do k = 1, size(s_list)
         coefficients(1, :, k) = surfaces(k)%rmnc
         coefficients(2, :, k) = surfaces(k)%rmns
         coefficients(3, :, k) = surfaces(k)%zmnc
         coefficients(4, :, k) = surfaces(k)%zmns
         coefficients(5, :, k) = surfaces(k)%vmnc
         coefficients(6, :, k) = surfaces(k)%vmns
         coefficients(7, :, k) = surfaces(k)%bmnc
         coefficients(8, :, k) = surfaces(k)%bmns
         fluxes(:, k) = [surfaces(k)%iota, surfaces(k)%curr_i, surfaces(k)%curr_g]
         physical(:, k) = [surfaces(k)%pp_file, surfaces(k)%sqrtg00]
      end do
      state = [equilibrium%r0, equilibrium%z0, equilibrium%psitor_max, &
         equilibrium%sig_theta]
      rminor = equilibrium%rminor
      if (check_finite) then
         if (any(.not. ieee_is_finite(coefficients)) .or. any(.not. ieee_is_finite(fluxes)) .or. &
             any(.not. ieee_is_finite(state)) .or. any(.not. ieee_is_finite(physical)) .or. &
             .not. ieee_is_finite(rminor)) then
            status_code = 7
            overridden = .false.
            message = 'EFIT contour transform produced nonfinite output'
            call python_error(trim(message))
            return
         end if
      end if
      if (.not. overridden) message = ''
   end subroutine contour_transform_worker

   subroutine init_legacy

      integer :: nstep, nsurfmax, nsurf, nt, mpol, inp_label
      double precision :: hs, htheta, oneovernt, twoovernt, sigma

      call deinit

      open (1, file='efit_to_boozer.inp', status='old', action='read')
      read (1, *) nstep    !number of integration steps
      read (1, *) nlabel   !grid size over radial variable
      read (1, *) ntheta   !grid size over poloidal angle
      read (1, *) nsurfmax !number of starting points between the
      !magnetic axis and right box boundary
      !when searching for the separatrix
      read (1, *) nsurf    !number of flux surfaces in Boozer file
      read (1, *) mpol     !number of poloidal modes in Boozer file
      read (1, *) psimax   !psi at plasma boundary
      close (1)

      allocate (rbeg(nlabel), rsmall(nlabel), qsaf(nlabel), psi_pol(0:nlabel))
      allocate (psi_tor_vac(nlabel), psi_tor(0:nlabel), C_const(nlabel))
      allocate (R_spl(0:nspl, 0:ntheta, nlabel), Z_spl(0:nspl, 0:ntheta, nlabel), bmod_spl(0:nspl, 0:ntheta, nlabel))
      allocate (sqgnorm_spl(0:nspl, 0:ntheta, nlabel), Gfunc_spl(0:nspl, 0:ntheta, nlabel))

      call flint_for_Boozer(nstep, nsurfmax, nlabel, ntheta, &
                            rmn, rmx, zmn, zmx, raxis, zaxis, sigma, &
                            rbeg, rsmall, qsaf, &
                            psi_pol(1:nlabel), psi_tor_vac, &
                            psi_tor(1:nlabel), C_const, &
                            R_spl(0, 1:ntheta, :), &
                            Z_spl(0, 1:ntheta, :), &
                            bmod_spl(0, 1:ntheta, :), &
                            sqgnorm_spl(0, 1:ntheta, :), &
                            Gfunc_spl(0, 1:ntheta, :))

      call spline_magdata_in_symfluxcoord

      print *, 'Splining done'

      nt = ntheta

      inp_label = 1
      hs = 1.d0/dfloat(nsurf)
      oneovernt = 1.d0/dfloat(nt)
      twoovernt = 2.d0*oneovernt
      htheta = twopi*oneovernt
      allocate (calE(nt), calEm_times_Jb(nt))
      allocate (R_oft(nt), Z_oft(nt), al_oft(nt), B_oft(nt))
      allocate (Rmn_c(0:mpol), Rmn_s(0:mpol), Zmn_c(0:mpol), Zmn_s(0:mpol))
      allocate (almn_c(0:mpol), almn_s(0:mpol), Bmn_c(0:mpol), Bmn_s(0:mpol))

   end subroutine init_legacy

   subroutine deinit

      initialized = .false.
      icall = 0
      icall_c = 0
      call reset_field_eq_state()

      if (allocated(rbeg)) deallocate (rbeg)
      if (allocated(rsmall)) deallocate (rsmall)
      if (allocated(qsaf)) deallocate (qsaf)
      if (allocated(psi_pol)) deallocate (psi_pol)
      if (allocated(psi_tor_vac)) deallocate (psi_tor_vac)
      if (allocated(psi_tor)) deallocate (psi_tor)
      if (allocated(C_const)) deallocate (C_const)
      if (allocated(R_spl)) deallocate (R_spl)
      if (allocated(Z_spl)) deallocate (Z_spl)
      if (allocated(bmod_spl)) deallocate (bmod_spl)
      if (allocated(sqgnorm_spl)) deallocate (sqgnorm_spl)
      if (allocated(Gfunc_spl)) deallocate (Gfunc_spl)
      if (allocated(calE)) deallocate (calE)
      if (allocated(calEm_times_Jb)) deallocate (calEm_times_Jb)
      if (allocated(R_oft)) deallocate (R_oft)
      if (allocated(Z_oft)) deallocate (Z_oft)
      if (allocated(al_oft)) deallocate (al_oft)
      if (allocated(B_oft)) deallocate (B_oft)
      if (allocated(Rmn_c)) deallocate (Rmn_c)
      if (allocated(Rmn_s)) deallocate (Rmn_s)
      if (allocated(Zmn_c)) deallocate (Zmn_c)
      if (allocated(Zmn_s)) deallocate (Zmn_s)
      if (allocated(almn_c)) deallocate (almn_c)
      if (allocated(almn_s)) deallocate (almn_s)
      if (allocated(Bmn_c)) deallocate (Bmn_c)
      if (allocated(Bmn_s)) deallocate (Bmn_s)

   end subroutine deinit

   subroutine magdata(inp_label, s, psi, theta, q, dq_ds, C_norm, dC_norm_ds, &
                      sqrtg, bmod, dbmod_dtheta, R, dR_ds, dR_dtheta, &
                      Z, dZ_ds, dZ_dtheta, G, dG_ds, dG_dtheta)

      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      use efit_to_boozer_mod, only: psipol_max
      integer, intent(in) :: inp_label
      double precision, intent(inout) :: s, psi, theta
      double precision, intent(out) :: q, dq_ds, C_norm, dC_norm_ds, sqrtg, bmod, &
         dbmod_dtheta, R, dR_ds, dR_dtheta, Z, dZ_ds, dZ_dtheta, G, dG_ds, dG_dtheta

      if (.not. initialized) then
         call python_error('EFIT is not initialized; call init() first')
         return
      end if
      if (inp_label /= 1 .and. inp_label /= 2) then
         call python_error('EFIT inp_label must be 1 or 2')
         return
      end if
      if (.not. ieee_is_finite(theta)) then
         call python_error('EFIT theta must be finite')
         return
      end if
      if (inp_label == 1) then
         if (.not. ieee_is_finite(s) .or. s < 0.0d0 .or. s > 1.0d0) then
            call python_error('EFIT normalized toroidal flux must lie in [0, 1]')
            return
         end if
      else
         if (.not. ieee_is_finite(psi)) then
            call python_error('EFIT poloidal flux must be finite')
            return
         end if
         if (psi/psipol_max < 0.0d0 .or. psi/psipol_max > 1.0d0) then
            call python_error('EFIT poloidal flux lies outside the initialized domain')
            return
         end if
      end if
      call magdata_in_symfluxcoord_ext(inp_label, s, psi, theta, q, dq_ds, C_norm, dC_norm_ds, &
                                       sqrtg, bmod, dbmod_dtheta, R, dR_ds, dR_dtheta, Z, dZ_ds, dZ_dtheta, G, dG_ds, dG_dtheta)
   end subroutine magdata

end module efit_to_boozer
