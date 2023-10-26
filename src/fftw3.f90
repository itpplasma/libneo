! wrap FFTW3 bindings in module
module FFTW3
  use, intrinsic :: iso_c_binding
  include 'fftw3.f03'
end module FFTW3
