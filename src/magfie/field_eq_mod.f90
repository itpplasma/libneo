module field_eq_mod
   implicit none

   integer, parameter :: dp = kind(1.0d0)

   logical :: use_fpol = .true.                                      !<=18.12.18
   logical :: skip_read = .false.
   integer :: icall_eq = 0
   integer :: nrad, nzet, icp, nwindow_r, nwindow_z
   integer :: ierrfield

   ! When .true., magfie accepts guiding-center points outside the separatrix
   ! (scrape-off layer) instead of flagging ierrfield=1. The equilibrium field
   ! is splined over the whole EQDSK box and F is frozen to its vacuum value
   ! there, so it stays valid. Default .false. keeps the historic domain bound
   ! for all other consumers; edge-orbit tracing (POTATO) sets it .true.
   logical :: allow_sol = .false.

   real(dp) :: psib, btf, rtf, hrad, hzet
   real(dp) :: psi_axis, psi_sep, hfpol                            !<=18.12.18
   real(dp), dimension(:, :), allocatable :: psi, psi0
   real(dp), dimension(:, :), allocatable :: splfpol               !<=18.12.18
   real(dp), dimension(:, :, :), allocatable :: splpsi
   real(dp), dimension(:), allocatable :: rad, zet, xi, f
   integer, dimension(:), allocatable :: imi, ima, jmi, jma
   integer, dimension(:, :), allocatable :: ipoint

   ! Field variables for POTATO integration
   real(dp) :: psif, dpsidr, dpsidz, d2psidr2, d2psidrdz, d2psidz2

   ! Make temporary variables threadprivate
   !$omp threadprivate(psif,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2)

contains

   subroutine reset_field_eq_state()
      implicit none

      if (allocated(psi)) deallocate (psi)
      if (allocated(psi0)) deallocate (psi0)
      if (allocated(splfpol)) deallocate (splfpol)
      if (allocated(splpsi)) deallocate (splpsi)
      if (allocated(rad)) deallocate (rad)
      if (allocated(zet)) deallocate (zet)
      if (allocated(xi)) deallocate (xi)
      if (allocated(f)) deallocate (f)
      if (allocated(imi)) deallocate (imi)
      if (allocated(ima)) deallocate (ima)
      if (allocated(jmi)) deallocate (jmi)
      if (allocated(jma)) deallocate (jma)
      if (allocated(ipoint)) deallocate (ipoint)

      skip_read = .false.
      icall_eq = 0
   end subroutine reset_field_eq_state

end module field_eq_mod
