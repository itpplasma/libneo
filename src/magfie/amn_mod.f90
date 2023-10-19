module amn_mod
  use libneo_kinds, only : complex_kind, real_kind

  implicit none

  ! Fourier ampitudes of the original field:
  integer :: ntor_amn=1,mpol,ntor_ff,mpol_ff,nsqpsi,icall=0
  real(kind=real_kind) :: sqpsimin,sqpsimax,hsqpsi
  complex(kind=complex_kind), dimension(:,:,:,:), allocatable :: splapsi,splatet
  complex(kind=complex_kind), dimension(:,:), allocatable :: amnpsi, amntet, &
                                                 amnpsi_s, amntet_s,   &
                                                 amnpsi_ss,amntet_ss
  complex(kind=complex_kind), dimension(:), allocatable :: expthe,expphi
  ! Formfactors:
  integer :: nsqpsi_ff,nmodes_ff
  real(kind=real_kind) :: sqpsimin_ff,sqpsimax_ff,hsqpsi_ff
  integer,        dimension(:,:), allocatable :: ipoi_ff
  complex(kind=complex_kind), dimension(:,:,:), allocatable :: splffp,splfft
  complex(kind=complex_kind), dimension(:), allocatable :: fmnpsi, fmntet, &
                                                 fmnpsi_s, fmntet_s,   &
                                                 fmnpsi_ss,fmntet_ss
end module amn_mod
