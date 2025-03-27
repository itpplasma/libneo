module amn_mod
  use libneo_kinds, only : cdp, dp

  implicit none

  ! Fourier ampitudes of the original field:
  integer :: ntor_amn=1,mpol,ntor_ff,mpol_ff,nsqpsi,icall=0
  real(dp) :: sqpsimin,sqpsimax,hsqpsi
  complex(cdp), dimension(:,:,:,:), allocatable :: splapsi,splatet
  complex(cdp), dimension(:,:), allocatable :: amnpsi, amntet, &
                                                 amnpsi_s, amntet_s,   &
                                                 amnpsi_ss,amntet_ss
  complex(cdp), dimension(:), allocatable :: expthe,expphi
  ! Formfactors:
  integer :: nsqpsi_ff,nmodes_ff
  real(dp) :: sqpsimin_ff,sqpsimax_ff,hsqpsi_ff
  integer,        dimension(:,:), allocatable :: ipoi_ff
  complex(cdp), dimension(:,:,:), allocatable :: splffp,splfft
  complex(cdp), dimension(:), allocatable :: fmnpsi, fmntet, &
                                                 fmnpsi_s, fmntet_s,   &
                                                 fmnpsi_ss,fmntet_ss
end module amn_mod
