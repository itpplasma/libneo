from os.path import expanduser
from fffi import FortranLibrary, FortranModule

libmc_efit = FortranLibrary('mc_efit', path=expanduser('~/src/libneo/build'))
field_eq = FortranModule(libmc_efit, name='field_eq_mod')

# TODO: add more definitions useful for testing
# Attention: private variables and routines are inaccessible

field_eq.fdef("""
  real(8) :: psif, psib
""")

libmc_efit.fdef("""
  subroutine field(r,p,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
  
  real(8), intent(in) :: r,p,z
  real(8), intent(out) :: Br,Bp,Bz,dBrdR,dBrdp,dBrdZ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
  
  end subroutine field
  
  subroutine magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
    double precision, dimension(3),       intent(in)         :: x
    double precision,                     intent(out)        :: bmod
    double precision,                     intent(out)        :: sqrtg
    double precision, dimension(3),       intent(out)        :: bder
    double precision, dimension(3),       intent(out)        :: hcovar
    double precision, dimension(3),       intent(out)        :: hctrvr
    double precision, dimension(3),       intent(out)        :: hcurl
  end

  subroutine velo(tau, z, vz)
    double precision              ,       intent(in)         :: tau
    double precision, dimension(3),       intent(in)         :: z
    double precision, dimension(3),       intent(out)        :: vz
  end
""")

libmc_efit.compile(verbose=1)
field_eq.load()

