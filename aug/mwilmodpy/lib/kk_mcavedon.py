import ctypes as ct
import numpy as np

# Input: byref for arrays and strings, ctypes for scalars
# Output: always byref

if ct.sizeof(ct.c_long) == 8:
  libso = 'lib64'
else:
  libso = 'lib'
libkk = ct.cdll.LoadLibrary('/afs/ipp/aug/ads/'+libso+'/@sys/libkk8.so')


exp_eq='AUGD'
dia_eq='EQH'

class kkhelp:
  status = False

class KK:

  def __init__(self):
    self.ed = 0

  def kkeqints(self,isw):

    # Input
    c_err=ct.c_long(0)
    _err=ct.byref(c_err)
    c_isw=ct.c_long(isw)
    _isw=ct.byref(c_isw)

    status=libkk.kkeqints(_err,c_isw)
    return status

  def const_ed(self,nshot,exp=exp_eq,diag=dia_eq):
    out = self.kkeqtp(nshot,exp=exp_eq,diag=dia_eq,ed=0)
    self.ed = out.ed
    print('Using '+exp_eq+':'+dia_eq+' edition '+str(self.ed))

  def kkeqtp(self,nshot,exp=exp_eq,diag=dia_eq,ed=0):
    if ed == 0:
        if hasattr(self,'ed'):
            ed = self.ed

    ntim=20000
    # Input
    c_ed=ct.c_long(ed)
    _ed=ct.byref(c_ed)
    c_err=ct.c_long(0)
    _err=ct.byref(c_err)
    c_shot=ct.c_long(nshot)
    _shot=ct.byref(c_shot)
    c_exp=ct.c_char_p(exp.encode('utf-8'))
    c_dia=ct.c_char_p(diag.encode('utf-8'))
    # Special input
    c_nprnt=ct.c_long(0)
    _nprnt=ct.byref(c_nprnt)
    # Output
    c_ntim=ct.c_int(ntim)
    _ntim=ct.byref(c_ntim)
    tim_kk=(ct.c_float * ntim)()
    _tim_kk=ct.byref(tim_kk)

    status=libkk.kkeqtp(_err,c_exp,c_dia,c_shot,_ed,_nprnt,_ntim,_tim_kk)

    output = kkhelp()
    output.ed  = c_ed.value
    output.err = c_err.value
    #  Numpy output
    nt=c_ntim.value
    output.time=np.zeros(nt)
    for jt in range(nt):
      output.time[jt]=tim_kk[jt]
    return output

  def kkrzbrzt(self,nshot,tshot,Rin,zin,exp=exp_eq,diag=dia_eq,ed=0):
    if ed == 0:
        if hasattr(self,'ed'):
            ed = self.ed

    # Input
    nr=len(Rin)
    c_nr=ct.c_long(nr)
    _nr=ct.byref(c_nr)
    c_ed=ct.c_long(ed)
    _ed=ct.byref(c_ed)
    c_err=ct.c_long(0)
    _err=ct.byref(c_err)
    c_shot=ct.c_long(nshot)
    _shot=ct.byref(c_shot)
    c_exp=ct.c_char_p(exp.encode('utf-8'))
    c_dia=ct.c_char_p(diag.encode('utf-8'))
    c_tshot=ct.c_float(tshot)
    _tshot=ct.byref(c_tshot)
    # Special input
    c_Rin=(ct.c_float * nr)()
    c_Rin[0:nr]=Rin[0:nr]
    _Rin=ct.byref(c_Rin)
    c_zin=(ct.c_float * nr)()
    c_zin[0:nr]=zin[0:nr]
    _zin=ct.byref(c_zin)
    # Output
    br  = (ct.c_float * nr)()
    bz  = (ct.c_float * nr)()
    bt  = (ct.c_float * nr)()
    fpf = (ct.c_float * nr)()
    fjf = (ct.c_float * nr)()
    _br = ct.byref(br)
    _bz = ct.byref(bz)
    _bt = ct.byref(bt)
    _fpf = ct.byref(fpf)
    _fjf = ct.byref(fjf)

    status=libkk.kkrzbrzt(_err,c_exp,c_dia,c_shot,_ed,_tshot,_Rin,_zin,c_nr, \
                           _br,_bz,_bt,_fpf,_fjf)

    output = kkhelp()
    output.ed  = c_ed.value
    output.err = c_err.value
    # Numpy output
    output.br=np.zeros(nr)
    output.bz=np.zeros(nr)
    output.bt=np.zeros(nr)
    output.pf=np.zeros(nr)
    output.jf=np.zeros(nr)
    for jr in range (nr):
      output.br[jr]=br[jr]
      output.bz[jr]=bz[jr]
      output.bt[jr]=bt[jr]
      output.pf[jr]=fpf[jr]
      output.jf[jr]=fjf[jr]

    return output

  def kkrzptfn(self,nshot,tshot,Rin,zin,exp=exp_eq,diag=dia_eq,ed=0):
    if ed == 0:
        if hasattr(self,'ed'):
            ed = self.ed

    # Input
    nr=len(Rin)
    c_nr=ct.c_long(nr)
    _nr=ct.byref(c_nr)
    c_ed=ct.c_long(ed)
    _ed=ct.byref(c_ed)
    c_err=ct.c_long(0)
    _err=ct.byref(c_err)
    c_shot=ct.c_long(nshot)
    _shot=ct.byref(c_shot)
    c_exp=ct.c_char_p(exp.encode('utf-8'))
    c_dia=ct.c_char_p(diag.encode('utf-8'))
    c_tshot=ct.c_float(tshot)
    _tshot=ct.byref(c_tshot)
    # Special input
    c_Rin=(ct.c_float * nr)()
    c_Rin[0:nr]=Rin[0:nr]
    _Rin=ct.byref(c_Rin)
    c_zin=(ct.c_float * nr)()
    c_zin[0:nr]=zin[0:nr]
    _zin=ct.byref(c_zin)
    # Output
    fpf = (ct.c_float * nr)()
    ftf = (ct.c_float * nr)()
    rhop= (ct.c_float * nr)()
    rhot= (ct.c_float * nr)()
    _fpf = ct.byref(fpf)
    _ftf = ct.byref(ftf)
    _rhop= ct.byref(rhop)
    _rhot= ct.byref(rhot)

    status=libkk.kkrzptfn(_err,c_exp,c_dia,c_shot,_ed,_tshot,_Rin,_zin,c_nr, \
                           _fpf,_rhop,_ftf,_rhot)

    output = kkhelp()
    output.ed  = c_ed.value
    output.err = c_err.value
    # Numpy output
    output.pf=np.zeros(nr)
    output.tf=np.zeros(nr)
    output.rho_p=np.zeros(nr)
    output.rho_t=np.zeros(nr)
    for jr in range (nr):
      output.pf[jr]=fpf[jr]
      output.tf[jr]=ftf[jr]
      output.rho_p[jr]=rhop[jr]
      output.rho_t[jr]=rhot[jr]

    return output

  def kkrhopto(self,nshot,tshot,rhop,exp=exp_eq,diag=dia_eq,ed=0):
    if ed == 0:
        if hasattr(self,'ed'):
            ed = self.ed

    nr=len(rhop)

    # Input
    c_nr=ct.c_long(nr)
    _nr=ct.byref(c_nr)
    c_ed=ct.c_long(ed)
    _ed=ct.byref(c_ed)
    c_err=ct.c_long(0)
    _err=ct.byref(c_err)
    c_shot=ct.c_long(nshot)
    _shot=ct.byref(c_shot)
    c_exp=ct.c_char_p(exp.encode('utf-8'))
    c_dia=ct.c_char_p(diag.encode('utf-8'))
    c_tshot=ct.c_float(tshot)
    _tshot=ct.byref(c_tshot)
    # Special input
    c_rhop=(ct.c_float * nr)()
    c_rhop[0:nr]=rhop[0:nr]
    _rhop=ct.byref(c_rhop)
    lexp = ct.c_long( len(exp) )
    ldia = ct.c_long( len(diag) )
    # Output
    fpf = (ct.c_float * nr)()
    ftf = (ct.c_float * nr)()
    rhot= (ct.c_float * nr)()
    _fpf = ct.byref(fpf)
    _ftf = ct.byref(ftf)
    _rhot= ct.byref(rhot)

    status=libkk.kkrhopto(_err,c_exp,c_dia,c_shot,_ed,_tshot,_rhop,c_nr, \
                          _rhot,_fpf,_ftf)

    output = kkhelp()
    output.ed  = c_ed.value
    output.err = c_err.value
    # Numpy output
    output.rho_t=np.zeros(nr)
    output.pf=np.zeros(nr)
    output.tf=np.zeros(nr)
    for jr in range(nr):
      output.rho_t[jr]= rhot[jr]
      output.pf[jr]     = fpf[jr]
      output.tf[jr]     = ftf[jr]

    return output

  def kkeqpsp(self,nshot,tshot,flux,exp=exp_eq,diag=dia_eq,ed=0):
    if ed == 0:
        if hasattr(self,'ed'):
            ed = self.ed

    # Input
    c_ed=ct.c_long(ed)
    c_err=ct.c_long(0)
    _err=ct.byref(c_err)
    _ed=ct.byref(c_ed)
    c_shot=ct.c_long(nshot)
    _shot=ct.byref(c_shot)
    c_exp=ct.c_char_p(exp.encode('utf-8'))
    c_dia=ct.c_char_p(diag.encode('utf-8'))
    c_tshot=ct.c_float(tshot)
    _tshot=ct.byref(c_tshot)
    c_tsfh=ct.c_float(0.)
    _tsfh=ct.byref(c_tsfh)
    lexp = ct.c_long( len(exp) )
    ldia = ct.c_long( len(diag) )
    # Specific input
    c_flux=ct.c_float(flux)
    _flux=ct.byref(c_flux)
    # Output
    nrzmax=1000
    npfm=1001
    mpfm=1001
    pfm=(ct.c_float * mpfm * npfm )()
    _pfm=ct.byref(pfm)
    c_ndim=ct.c_long(1000)
    _ndim=ct.byref(c_ndim)
    rz=(ct.c_float * nrzmax * 2)()
    _rz=ct.byref(rz)
    c_nrz=ct.c_long(nrzmax)
    _nrz=ct.byref(c_nrz)
    status=libkk.kkeqpsp(_err,c_exp,c_dia,c_shot,_ed,_tshot, \
                         c_flux,_ndim,_rz,_nrz,_pfm,_tsfh)

    output = kkhelp()
    output.ed  = c_ed.value
    output.err = c_err.value
    # Numpy output
    n_rz=c_nrz.value
    if (n_rz > 0):
      output.r_surf=np.zeros(n_rz)
      output.z_surf=np.zeros(n_rz)
      for jrz in range(n_rz):
        output.r_surf[jrz]=rz[0][jrz]
        output.z_surf[jrz]=rz[1][jrz]
    else:
      output.r_surf=np.zeros(1)
      output.z_surf=np.zeros(1)
    return output

  def kkeqtfl(self,nshot,tshot,exp=exp_eq,diag=dia_eq,ed=0):
    if ed == 0:
        if hasattr(self,'ed'):
            ed = self.ed

    ndxymax = 3200
    # Input
    c_ed=ct.c_long(ed)
    _ed=ct.byref(c_ed)
    c_err=ct.c_long(0)
    _err=ct.byref(c_err)
    c_shot=ct.c_long(nshot)
    _shot=ct.byref(c_shot)
    c_exp=ct.c_char_p(exp.encode('utf-8'))
    c_dia=ct.c_char_p(diag.encode('utf-8'))
    c_tshot=ct.c_float(tshot)
    _tshot=ct.byref(c_tshot)
    # Output
    c_pfl = (ct.c_float * ndxymax)()
    c_tfl = (ct.c_float * ndxymax)()
    c_lpf = ct.c_long(ndxymax)
    _pfl  = ct.byref(c_pfl)
    _tfl  = ct.byref(c_tfl)
    _lpf  = ct.byref(c_lpf)

    status = libkk.kkeqtfl(_err,c_exp,c_dia,c_shot,_ed,_tshot, \
                           _lpf,_pfl,_tfl)

    output = kkhelp()
    output.ed  = c_ed.value
    output.err = c_err.value
    # Numpy output
    if c_err.value !=0:
      print(' ERROR: kkeqtfl')
    else:
      nr = c_lpf.value+1
      output.pf = np.zeros(nr,dtype=float)
      output.tf = np.zeros(nr,dtype=float)
      for jr in range(nr):
        output.pf[jr] = c_pfl[jr]
        output.tf[jr] = c_tfl[jr]
    return output

  def kkpfrhop(self,nshot,tshot,fpf,exp=exp_eq,diag=dia_eq,ed=0):
    if ed == 0:
        if hasattr(self,'ed'):
            ed = self.ed

    nr = len(fpf)
    # Input
    c_ed=ct.c_long(ed)
    _ed=ct.byref(c_ed)
    c_err=ct.c_long(0)
    _err=ct.byref(c_err)
    c_shot=ct.c_long(nshot)
    _shot=ct.byref(c_shot)
    c_exp=ct.c_char_p(exp.encode('utf-8'))
    c_dia=ct.c_char_p(diag.encode('utf-8'))
    c_tshot=ct.c_float(tshot)
    _tshot=ct.byref(c_tshot)
    # Special input
    c_fp = (ct.c_float * nr)()
    c_fp[0:nr]=fpf[0:nr]
    _fp = ct.byref(c_fp)
    c_nr = ct.c_long(nr)
    _nr=ct.byref(c_nr)
    # Output
    c_rhop = (ct.c_float * nr)()
    _rhop = ct.byref(c_rhop)

    status=libkk.kkpfrhop(_err,c_exp,c_dia,c_shot,_ed,_tshot, \
                          _fp,c_nr,_rhop)

    output = kkhelp()
    output.ed  = c_ed.value
    output.err = c_err.value
    # Numpy output
    output.rho_p=np.zeros(nr)
    for jr in range(nr):
      output.rho_p[jr]=c_rhop[jr]
    return output

  def kktfrhot(self,nshot,tshot,ftf,exp=exp_eq,diag=dia_eq,ed=0):
    if ed == 0:
        if hasattr(self,'ed'):
            ed = self.ed

    nr = len(ftf)
    # Input
    c_ed=ct.c_long(ed)
    _ed=ct.byref(c_ed)
    c_err=ct.c_long(0)
    _err=ct.byref(c_err)
    c_shot=ct.c_long(nshot)
    _shot=ct.byref(c_shot)
    c_exp=ct.c_char_p(exp.encode('utf-8'))
    c_dia=ct.c_char_p(diag.encode('utf-8'))
    c_tshot=ct.c_float(tshot)
    _tshot=ct.byref(c_tshot)
    # Special input
    c_ft = (ct.c_float * nr)()
    c_ft[0:nr]=ftf[0:nr]
    _ft = ct.byref(c_ft)
    c_nr = ct.c_long(nr)
    _nr=ct.byref(c_nr)
    # Output
    c_rhot = (ct.c_float * nr)()
    _rhot = ct.byref(c_rhot)

    status=libkk.kktfrhot(_err,c_exp,c_dia,c_shot,_ed,_tshot, \
                          _ft,c_nr,_rhot)

    output = kkhelp()
    output.ed  = c_ed.value
    output.err = c_err.value
    # Numpy output
    output.rho_t=np.zeros(nr)
    for jr in range(nr):
      output.rho_t[jr]=c_rhot[jr]
    return output

  def kkrhopfa(self,nshot,tshot,rhop,exp=exp_eq,diag=dia_eq,ed=0):
    if ed == 0:
        if hasattr(self,'ed'):
            ed = self.ed

    nr=len(rhop)
    # Input
    c_ed=ct.c_long(ed)
    _ed=ct.byref(c_ed)
    c_err=ct.c_long(0)
    _err=ct.byref(c_err)
    c_shot=ct.c_long(nshot)
    _shot=ct.byref(c_shot)
    c_exp=ct.c_char_p(exp.encode('utf-8'))
    c_dia=ct.c_char_p(diag.encode('utf-8'))
    c_tshot=ct.c_float(tshot)
    _tshot=ct.byref(c_tshot)
    # Special input
    c_nr = ct.c_long(nr)
    _nr = ct.byref(c_nr)
    c_rhop = (ct.c_float * nr)()
    c_rhop[0:nr]=rhop[0:nr]
    _rhop = ct.byref(c_rhop)
    lexp = ct.c_long( len(exp) )
    ldia = ct.c_long( len(diag) )
    # Output
    c_farea =  (ct.c_float * nr)()
    c_fpf =  (ct.c_float * nr)()
    _farea = ct.byref(c_farea)
    _fpf = ct.byref(c_fpf)

    status=libkk.kkrhopfa(_err,c_exp,c_dia,c_shot,_ed,_tshot, \
                          _rhop,c_nr,_farea,_fpf)

    output = kkhelp()
    output.ed  = c_ed.value
    output.err = c_err.value
    # Numpy output
    output.area=np.zeros(nr)
    output.pf=np.zeros(nr)
    for jr in range(nr):
      output.area[jr]=c_farea[jr]
      output.pf[jr]=c_fpf[jr]
    return output

  def kkrhopcs(self,nshot,tshot,rhop,exp=exp_eq,diag=dia_eq,ed=0):
    if ed == 0:
        if hasattr(self,'ed'):
            ed = self.ed

    nr=len(rhop)
    # Input
    c_ed=ct.c_long(ed)
    _ed=ct.byref(c_ed)
    c_err=ct.c_long(0)
    _err=ct.byref(c_err)
    c_shot=ct.c_long(nshot)
    _shot=ct.byref(c_shot)
    c_exp=ct.c_char_p(exp.encode('utf-8'))
    c_dia=ct.c_char_p(diag.encode('utf-8'))
    c_tshot=ct.c_float(tshot)
    _tshot=ct.byref(c_tshot)
    # Special input
    c_nr = ct.c_long(nr)
    _nr = ct.byref(c_nr)
    c_rhop = (ct.c_float * nr)()
    c_rhop[0:nr]=rhop[0:nr]
    _rhop = ct.byref(c_rhop)
    lexp = ct.c_long( len(exp) )
    ldia = ct.c_long( len(diag) )
    # Output
    c_fcirc =  (ct.c_float * nr)()
    c_fsurf =  (ct.c_float * nr)()
    c_fpf =  (ct.c_float * nr)()
    _fcirc = ct.byref(c_fcirc)
    _fsurf = ct.byref(c_fsurf)
    _fpf = ct.byref(c_fpf)

    status=libkk.kkrhopcs(_err,c_exp,c_dia,c_shot,_ed,_tshot, \
                          _rhop,c_nr,_fcirc,_fsurf,_fpf)

    output = kkhelp()
    output.ed  = c_ed.value
    output.err = c_err.value
    # Numpy output
    output.surf=np.zeros(nr)
    output.circ=np.zeros(nr)
    output.pf=np.zeros(nr)
    for jr in range(nr):
      output.circ[jr]=c_fcirc[jr]
      output.surf[jr]=c_fsurf[jr]
      output.pf[jr]=c_fpf[jr]
    return output

  def kkrhotpq(self,nshot,tshot,rhot,exp=exp_eq,diag=dia_eq,ed=0):
    if ed == 0:
        if hasattr(self,'ed'):
            ed = self.ed

    nr=len(rhot)
    # Input
    c_ed=ct.c_long(ed)
    _ed=ct.byref(c_ed)
    c_err=ct.c_long(0)
    _err=ct.byref(c_err)
    c_shot=ct.c_long(nshot)
    _shot=ct.byref(c_shot)
    c_exp=ct.c_char_p(exp.encode('utf-8'))
    c_dia=ct.c_char_p(diag.encode('utf-8'))
    c_tshot=ct.c_float(tshot)
    _tshot=ct.byref(c_tshot)
    # Special input
    c_nr = ct.c_long(nr)
    _nr = ct.byref(c_nr)
    c_rhot = (ct.c_float * nr)()
    c_rhot[0:nr]=rhot[0:nr]
    _rhot = ct.byref(c_rhot)
    lexp = ct.c_long( len(exp) )
    ldia = ct.c_long( len(diag) )
    # Output
    c_qp =  (ct.c_float * nr)()
    c_fp =  (ct.c_float * nr)()
    c_ft =  (ct.c_float * nr)()
    _qp = ct.byref(c_qp)
    _fp = ct.byref(c_fp)
    _ft = ct.byref(c_ft)

    status=libkk.kkrhotpq(_err,c_exp,c_dia,c_shot,_ed,_tshot, \
                          _rhot,c_nr,_qp,_fp,_ft)

    output = kkhelp()
    output.ed  = c_ed.value
    output.err = c_err.value
    # Numpy output
    output.fp=np.zeros(nr)
    output.ft=np.zeros(nr)
    output.qp=np.zeros(nr)
    for jr in range(nr):
      output.fp[jr]=c_fp[jr]
      output.ft[jr]=c_ft[jr]
      output.qp[jr]=c_qp[jr]
    return output

  def kkeqpfx(self,nshot,tshot,exp=exp_eq,diag=dia_eq,ed=0):
    if ed == 0:
        if hasattr(self,'ed'):
            ed = self.ed

    nr=4
    nr1=nr+1
    # Input
    c_ed=ct.c_long(ed)
    _ed=ct.byref(c_ed)
    c_err=ct.c_long(0)
    _err=ct.byref(c_err)
    c_shot=ct.c_long(nshot)
    _shot=ct.byref(c_shot)
    c_exp=ct.c_char_p(exp.encode('utf-8'))
    c_dia=ct.c_char_p(diag.encode('utf-8'))
    c_tshot=ct.c_float(tshot)
    _tshot=ct.byref(c_tshot)
    # Special input
    c_lpfx=ct.c_long(nr)
    _lpfx=ct.byref(c_lpfx)
    # Output
    c_pfx =(ct.c_float * nr1)()
    c_rpfx=(ct.c_float * nr1)()
    c_zpfx=(ct.c_float * nr1)()
    _pfx =ct.byref(c_pfx)
    _rpfx=ct.byref(c_rpfx)
    _zpfx=ct.byref(c_zpfx)

    status=libkk.kkeqpfx(_err,c_exp,c_dia,c_shot,_ed,_tshot, \
                         _lpfx,_pfx,_rpfx,_zpfx)

    output = kkhelp()
    output.ed  = c_ed.value
    output.err = c_err.value
    # Float output
    output.Raxis=c_rpfx[0]
    output.zaxis=c_zpfx[0]
    output.Rspx=c_rpfx[1]
    output.zspx=c_zpfx[1]
    output.Rlim=c_rpfx[2]
    output.zlim=c_zpfx[2]
    return output

  def kkeqpres(self,nshot,tshot,exp=exp_eq,diag=dia_eq,ed=0):
    if ed == 0:
        if hasattr(self,'ed'):
            ed = self.ed

    npres_max=100000
    # Input
    c_ed=ct.c_long(ed)
    _ed=ct.byref(c_ed)
    c_err=ct.c_long(0)
    _err=ct.byref(c_err)
    c_shot=ct.c_long(nshot)
    _shot=ct.byref(c_shot)
    c_exp=ct.c_char_p(exp.encode('utf-8'))
    c_dia=ct.c_char_p(diag.encode('utf-8'))
    c_tshot=ct.c_float(tshot)
    _tshot=ct.byref(c_tshot)
    # Output
    c_npres=ct.c_long(npres_max)
    _npres=ct.byref(c_npres)
    c_fp   =(ct.c_float * npres_max)()
    c_pres =(ct.c_float * npres_max)()
    c_presp=(ct.c_float * npres_max)()
    _fp   =ct.byref(c_fp)
    _pres =ct.byref(c_pres)
    _presp=ct.byref(c_presp)

    status=libkk.kkeqpres(_err,c_exp,c_dia,c_shot,_ed,_tshot, \
                          _npres,_fp,_pres,_presp)

    n_pres=c_npres.value+1
    output = kkhelp()
    output.ed  = c_ed.value
    output.err = c_err.value
    # Numpy output
    output.pres = np.zeros(n_pres)
    output.fp   = np.zeros(n_pres)
    for jpres in range(n_pres):
      output.pres[jpres] = c_pres[jpres]
      output.fp[jpres]   = c_fp[jpres]
    return output

  def kkeqqpl(self,nshot,tshot,exp=exp_eq,diag=dia_eq,ed=0):

    lpf_max=100000
    # Input
    c_ed=ct.c_long(ed)
    _ed=ct.byref(c_ed)
    c_err=ct.c_long(0)
    _err=ct.byref(c_err)
    c_shot=ct.c_long(nshot)
    _shot=ct.byref(c_shot)
    c_exp=ct.c_char_p(exp.encode('utf-8'))
    c_dia=ct.c_char_p(diag.encode('utf-8'))
    c_tshot=ct.c_float(tshot)
    _tshot=ct.byref(c_tshot)
    # Output
    c_lpf=ct.c_long(lpf_max)
    _lpf=ct.byref(c_lpf)
    c_pfl=(ct.c_float * lpf_max)()
    c_qpl=(ct.c_float * lpf_max)()
    _pfl=ct.byref(c_pfl)
    _qpl=ct.byref(c_qpl)

    status = libkk.kkeqqpl(_err,c_exp,c_dia,c_shot,_ed,_tshot, \
                           _lpf,_pfl,_qpl)

    output = kkhelp()
    output.ed  = c_ed.value
    output.err = c_err.value
    #  Numpy output
    nrho=c_lpf.value+1
    output.pfl=np.zeros(nrho)
    output.qpl=np.zeros(nrho)
    for jrho in range(nrho):
      output.pfl[jrho]=c_pfl[jrho]
      output.qpl[jrho]=c_qpl[jrho]
    return output

  def kkeqffs(self,nshot,tshot,exp=exp_eq,diag=dia_eq,ed=0):
    if ed == 0:
        if hasattr(self,'ed'):
            ed = self.ed

    nr_max = 100000
    # Input
    c_ed=ct.c_long(ed)
    _ed=ct.byref(c_ed)
    c_err=ct.c_long(0)
    _err=ct.byref(c_err)
    c_shot=ct.c_long(nshot)
    _shot=ct.byref(c_shot)
    c_exp=ct.c_char_p(exp.encode('utf-8'))
    c_dia=ct.c_char_p(diag.encode('utf-8'))
    c_tshot=ct.c_float(tshot)
    _tshot=ct.byref(c_tshot)
    # Special input
    c_typ = ct.c_long(11)
    _typ = ct.byref(c_typ)
    # Output
    c_nr=ct.c_long(nr_max)
    _nr = ct.byref(c_nr)
    c_fp     = (ct.c_float * nr_max)()
    c_rbphi  = (ct.c_float * nr_max)()
    c_drbphi = (ct.c_float * nr_max)()
    _fp    = ct.byref(c_fp)
    _rbphi = ct.byref(c_rbphi)
    _drbphi = ct.byref(c_drbphi)

    status=libkk.kkeqffs(_err,c_exp,c_dia,c_shot,_ed,_tshot, \
                         c_typ,_nr,_fp,_rbphi,_drbphi)

    n_rbphi=c_nr.value+1
    output = kkhelp()
    output.ed  = c_ed.value
    output.err = c_err.value
    # Numpy output
    output.fp   =np.zeros(n_rbphi)
    output.rb=np.zeros(n_rbphi)
    for jr in range(n_rbphi):
      output.fp[jr]=c_fp[jr]
      output.rb[jr]=c_rbphi[jr]
    return output

  def kkcutsx(self,nshot,tshot,Rin,zin,phi,exp=exp_eq,diag=dia_eq,ed=0):
    if ed == 0:
        if hasattr(self,'ed'):
            ed = self.ed

    # Input
    c_ed=ct.c_long(ed)
    _ed=ct.byref(c_ed)
    c_err=ct.c_long(0)
    _err=ct.byref(c_err)
    c_shot=ct.c_long(nshot)
    _shot=ct.byref(c_shot)
    c_exp=ct.c_char_p(exp.encode('utf-8'))
    c_dia=ct.c_char_p(diag.encode('utf-8'))
    c_tshot=ct.c_float(tshot)
    _tshot=ct.byref(c_tshot)
    # Special input
    c_Rin=ct.c_float(Rin)
    c_zin=ct.c_float(zin)
    c_phi=ct.c_float(phi)
    _Rin=ct.byref(c_Rin)
    _zin=ct.byref(c_zin)
    _phi=ct.byref(c_phi)
    #Output
    c_Rout =ct.c_float(0)
    c_zout =ct.c_float(0)
    c_psisx=ct.c_float(0)
    c_tsfh =ct.c_float(0)
    _Rout =ct.byref(c_Rout)
    _zout =ct.byref(c_zout)
    _psisx=ct.byref(c_psisx)
    _tsfh =ct.byref(c_tsfh)

    status=libkk.kkcutsx(_err,c_exp,c_dia,c_shot,_ed,_tshot, \
                         c_Rin,c_zin,c_phi,_Rout,_zout,_psisx,_tsfh)

    output = kkhelp()
    output.ed  = c_ed.value
    output.err = c_err.value
    output.Rsx = c_Rout.value
    output.zsx = c_zout.value
    return output

  def kkrhorz(self,nshot,tshot,rhop,angle=0,exp=exp_eq,diag=dia_eq,ed=0):
    if ed == 0:
        if hasattr(self,'ed'):
            ed = self.ed

    nr=len(rhop)

    # Input
    c_nr=ct.c_long(nr)
    _nr=ct.byref(c_nr)
    c_ed=ct.c_long(ed)
    _ed=ct.byref(c_ed)
    c_err=ct.c_long(0)
    _err=ct.byref(c_err)
    c_shot=ct.c_long(nshot)
    _shot=ct.byref(c_shot)
    c_exp=ct.c_char_p(exp.encode('utf-8'))
    c_dia=ct.c_char_p(diag.encode('utf-8'))
    c_tshot=ct.c_float(tshot)
    _tshot=ct.byref(c_tshot)
    # Special input
    c_iorg=ct.c_long(0)
    _iorg=ct.byref(c_iorg)
    c_rhop=(ct.c_float * nr)()
    c_rhop[0:nr]=rhop[0:nr]
    _rhop=ct.byref(c_rhop)
    lexp = ct.c_long( len(exp) )
    ldia = ct.c_long( len(diag) )
    c_angle=ct.c_float(angle)
    _angle=ct.byref(c_angle)
    # Output
    c_zn = (ct.c_float * nr)()
    c_rn = (ct.c_float * nr)()
    _zn = ct.byref(c_zn)
    _rn = ct.byref(c_rn)

    status=libkk.kkrhorz(_err,c_exp,c_dia,c_shot,_ed,_tshot, \
                         _rhop,c_nr,c_iorg,c_angle,_rn,_zn)

    output = kkhelp()
    output.ed  = c_ed.value
    output.err = c_err.value
    # Numpy output
    output.z=np.zeros(nr)
    output.r=np.zeros(nr)
    for jr in range(nr):
      output.z[jr]=c_zn[jr]
      output.r[jr]=c_rn[jr]

    return output

  def kkrhorz(self,nshot,tshot,rhop,angle=0,exp=exp_eq,diag=dia_eq,ed=0):
    if ed == 0:
        if hasattr(self,'ed'):
            ed = self.ed

    nr=len(rhop)

    # Input
    c_nr=ct.c_long(nr)
    _nr=ct.byref(c_nr)
    c_ed=ct.c_long(ed)
    _ed=ct.byref(c_ed)
    c_err=ct.c_long(0)
    _err=ct.byref(c_err)
    c_shot=ct.c_long(nshot)
    _shot=ct.byref(c_shot)
    c_exp=ct.c_char_p(exp.encode('utf-8'))
    c_dia=ct.c_char_p(diag.encode('utf-8'))
    c_tshot=ct.c_float(tshot)
    _tshot=ct.byref(c_tshot)
    # Special input
    c_iorg=ct.c_long(0)
    _iorg=ct.byref(c_iorg)
    c_rhop=(ct.c_float * nr)()
    c_rhop[0:nr]=rhop[0:nr]
    _rhop=ct.byref(c_rhop)
    lexp = ct.c_long( len(exp) )
    ldia = ct.c_long( len(diag) )
    c_angle=ct.c_float(angle)
    _angle=ct.byref(c_angle)
    # Output
    c_zn = (ct.c_float * nr)()
    c_rn = (ct.c_float * nr)()
    _zn = ct.byref(c_zn)
    _rn = ct.byref(c_rn)

    status=libkk.kkrhorz(_err,c_exp,c_dia,c_shot,_ed,_tshot, \
                         _rhop,c_nr,c_iorg,c_angle,_rn,_zn)

    output = kkhelp()
    output.ed  = c_ed.value
    output.err = c_err.value
    # Numpy output
    output.z=np.zeros(nr)
    output.r=np.zeros(nr)
    for jr in range(nr):
      output.z[jr]=c_zn[jr]
      output.r[jr]=c_rn[jr]

    return output

  def kkeqtpsa(self,nshot,tshot,divia='a',method=1,exp=exp_eq,diag=dia_eq,ed=0):
    """
    {FPP,EQU} calc. strike-points & dist, angle:
                kkEQTPsa (iERR  ,expnam,dianam,nSHOT,nEDIT, tSHOT,
               >                 DIVia,method,
               <                 RcTP,zcTP,scTP,acTP,       tSHf)
        DIVia  ...:= different target_plates etc. (char.*8)      ->
                  := 'TPLB' ...inner divertor Bottom (Div.I) 
                  := 'TPRB' ...outer divertor Bottom (Div.I)
                  := 'i'    ...inner divertor Bottom (Div.II & IIb)   
                  := 'a'    ...outer divertor Bottom (Div.II & IIb)
                  := 'ri'   ...roof baffle left part (Div.II & IIb)
                  := 'ra'   ...roof baffle right part(Div.II & IIb)
                  := 'TPLT' ...target_plate  Left_Top
                  := 'TPRT' ...target_plate Right_Top
        method ...method_#  at present == 1       (integer)      ->

        RcTP   ...R-coord. of strikepoint   [m]   (real *4)         <--
        zcTP   ...z-coord. of strikepoint   [m]   (real *4)         <--
        scTP   ...dist. of strikept [m] TP(Div.II)(real *4)         <--
        acTP   ...normal angle at strikepoint     (real *4) 
    """
    # Input
    c_ed=ct.c_long(ed)
    _ed=ct.byref(c_ed)
    c_err=ct.c_long(0)
    _err=ct.byref(c_err)
    c_shot=ct.c_long(nshot)
    _shot=ct.byref(c_shot)
    c_exp=ct.c_char_p(exp.encode('utf-8'))
    c_dia=ct.c_char_p(diag.encode('utf-8'))
    c_tshot=ct.c_float(tshot)
    _tshot=ct.byref(c_tshot)
    # Special input
    c_divia = ct.c_char_p(divia)
    c_method=ct.c_long(method)

    # Output
    nr = 1
    c_rctp = (ct.c_float * nr)()
    c_zctp = (ct.c_float * nr)()
    c_sctp = (ct.c_float * nr)()
    c_actp = (ct.c_float * nr)()
    c_tshf = (ct.c_float * nr)()
    _rctp = ct.byref(c_rctp)
    _zctp = ct.byref(c_zctp)
    _sctp = ct.byref(c_sctp)
    _actp = ct.byref(c_actp)
    _tshf = ct.byref(c_tshf)

    status = libkk.kkeqtpsa(_err, c_exp, c_dia, c_shot, _ed, _tshot,\
            c_divia, c_method,\
            _rctp, _zctp, _sctp, _actp, _tshf)

    output = kkhelp()
    output.ed  = c_ed.value
    output.err = c_err.value
    output.RcTP  = c_rctp[0]
    output.zcTP = c_zctp[0]
    output.scTP = c_sctp[0]
    output.acTP = c_actp[0]
    output.tSHf = c_tshf[0]

    return output

  def kkgcdimensions(self,nshot,exp=exp_eq,diag=dia_eq,ed=0):

    # Input
    c_ed=ct.c_long(ed)
    _ed=ct.byref(c_ed)
    c_err=ct.c_long(0)
    _err=ct.byref(c_err)
    c_shot=ct.c_long(nshot)
    _shot=ct.byref(c_shot)
    c_exp=ct.c_char_p(exp.encode('utf-8'))
    c_dia=ct.c_char_p(diag.encode('utf-8'))

    # Output
    c_shotgc = ct.c_long(0)
    _shotgc = ct.byref(c_shotgc)
    c_ndgc = ct.c_long(0)
    _ndgc = ct.byref(c_ndgc)
    c_ndit = ct.c_long(0)
    _ndit = ct.byref(c_ndit)
    c_ndim = ct.c_long(0)
    _ndim = ct.byref(c_ndim)
    c_nvalid = ct.c_long(0)
    _nvalid = ct.byref(c_nvalid)

    status = libkk.kkgcdimensions(_err, c_exp, c_dia, c_shot, _ed, \
            _shotgc, _ndgc, _ndit, _ndim, _nvalid)

    output = kkhelp()
    output.ed  = c_ed.value
    output.err = c_err.value
    output.shotgc = c_shotgc.value
    output.ndgc = c_ndgc.value
    output.ndit = c_ndit.value
    output.ndim = c_ndim.value
    output.nvalid = c_nvalid.value

    return output

  def kkgcread(self,nshot):
        exp = 'AUGD'
        diag = 'EQH'
        shot = nshot
        edition =ct.c_int(0)
        byref = ct.byref
        error = ct.c_int(0)
        shotgc = ct.c_int(0) # shot number of the corresponding YGC file
        ndgc = ct.c_int(0) # number of all structures (including unused ones)
        ndit = ct.c_int(0) # maximal number of points per single structure
        ndim = ct.c_int(0) # total number of points in all structures (including unused ones)
        nvalid = ct.c_int(0) # number of valied structures

        libkk.kkgcdimensions(byref(error), exp, diag, shot, byref(edition),
                     byref(shotgc), byref(ndgc), byref(ndit),
                     byref(ndim), byref(nvalid) )


        lenix = (ct.c_int*nvalid.value)()


        xygcs = (((ct.c_float*nvalid.value)*ndit.value)*2)() #numpy.zeros((ndgc.value,ndit.value))
        ngc = ct.c_int(0)
        libkk.kkgcread(ct.byref(error), exp, diag, shot, byref(edition),
               nvalid, ndit, ct.byref(xygcs), ct.byref(ngc), ct.byref(lenix) )

        rz = np.frombuffer(xygcs, np.float32).reshape(2, ndit.value, nvalid.value)
        lenix = np.frombuffer(lenix, np.int32)

        return [{'R':rz[0,:lenix[i],i], 'z':rz[1,:lenix[i],i]} for i in range(nvalid.value)]

  def fluxSurface(self,nshot,tshot,rhop=np.linspace(0.8,1.2,10),\
                  exp=exp_eq,diag=dia_eq,ed=0):

    rhop = np.append(rhop,1.)
    rhop = np.sort(rhop)
    output = self.kkrhopto(nshot,tshot,rhop,exp=exp_eq,diag=dia_eq,ed=ed)
    fluxSurface_r = []
    fluxSurface_z = []
    for flux in output.pf:
        output = self.kkeqpsp(nshot,tshot,flux,\
                              exp=exp_eq,diag=dia_eq,ed=ed)
        fluxSurface_r.append(output.r_surf)
        fluxSurface_z.append(output.z_surf)
    return fluxSurface_r,fluxSurface_z

  def fastRz2rho(self,nshot,R,z,timebase,exp=exp_eq,diag=dia_eq,ed=0,fast_dia='FPG',\
                 fast_exp='AUGD',fast_ed=0.,inORout='out',rho_pol=None):
    if nshot > 30150 and nshot < 31267 and diag=='EQH':
        exp = 'AUGE'
        fast_exp = 'AUGE'
    #read timebase of equilibrium
    out = self.kkeqtp(nshot,exp=exp,diag=diag,ed=ed)
    self.ed = out.ed
    eq_time = out.time
    #find time window between equilibrium and required mapping
    timebase = np.atleast_1d(timebase)
    jtBeg_eq = np.min(np.where(eq_time>timebase[0]))
    jtEnd_eq = np.max(np.where(eq_time<timebase[-1]))
    if timebase.size == 1:
        jtEnd_eq += 1
        jtBeg_eq -= 1
    #read raus/rin from fast diag
    from ddPlus import ddPlus as dd
    dd = dd.shotfile(fast_dia,nshot,experiment=fast_exp)
    if inORout == 'out':
        sig = dd('Raus')
        angle = 0.
    elif inORout == 'in':
        sig = dd('Rin')
        angle = 180.
    else:
        raise ValueError('inORout can only be in or out!')
    raus = sig.data
    raus_time = sig.time

    #calculate rho_pol and separatrix position (angle=0) for the time points of 
    #eq in between the required time mapping
    if_rho_pol = True
    if rho_pol is None:
        rho_pol = np.zeros((jtEnd_eq-jtBeg_eq,len(R)))
        if_rho_pol = False
    else:
        rho_pol = rho_pol[jtBeg_eq:jtEnd_eq]
    r_sep = np.zeros(jtEnd_eq-jtBeg_eq)
    r_sep_diff = np.array([])
    drho_dr = np.zeros(jtEnd_eq-jtBeg_eq)

    jtabs = 0
    for jt in range(jtBeg_eq,jtEnd_eq):
         if not if_rho_pol:
             rho_pol[jtabs] = self.kkrzptfn(nshot,eq_time[jt],R,z,exp=exp,diag=diag,ed=ed).rho_p
         out = self.kkrhorz(nshot,eq_time[jt],[0.99,1.],angle=angle,exp=exp,diag=diag,ed=ed)
         drho_dr[jtabs] = 0.01/np.diff(out.r)
         r_sep[jtabs] = out.r[1]
         #linear interpolate raus in between two time points of time_eq and calculate the difference with 
         #original signal. This would be the correction to the separatrix position of EQH
         #linar interp
         jtBeg_lin = np.argmin(np.abs(eq_time[jt]-raus_time))
         jtEnd_lin = np.argmin(np.abs(eq_time[jt+1]-raus_time))
         diff = np.interp(raus_time[jtBeg_lin:jtEnd_lin],\
                          [raus_time[jtBeg_lin],raus_time[jtEnd_lin]],\
                          [raus[jtBeg_lin],raus[jtEnd_lin]])-raus[jtBeg_lin:jtEnd_lin]
         r_sep_diff = np.append(diff,r_sep_diff)
         jtabs += 1

    jtBeg_lin = np.argmin(np.abs(eq_time[jtBeg_eq]-raus_time))
    jtEnd_lin = np.argmin(np.abs(eq_time[jtEnd_eq]-raus_time))
    #reinterpolate everything into requested timebase
    r_sep_diff = np.interp(timebase,raus_time[jtBeg_lin:jtEnd_lin],r_sep_diff)
    drho_dr = np.interp(timebase,eq_time[jtBeg_eq:jtEnd_eq],drho_dr)
    rho_corr = r_sep_diff*drho_dr
    rho_pol_end = np.zeros((len(timebase),rho_pol.shape[1]))
    for jr in range(rho_pol.shape[1]):
        rho_pol_end[:,jr] = np.interp(timebase,eq_time[jtBeg_eq:jtEnd_eq],rho_pol[:,jr])
        rho_pol_end[:,jr] += rho_corr
    return rho_pol_end

def main():
    kk = KK()
    import dd 
    sf = dd.shotfile('CNZ',31052,edition=1)
    R = sf.getAreaBase('R').data
    R = R[R != 0]
    z = sf.getAreaBase('z').data
    z = z[z != 0]
    phi = sf.getAreaBase('phi').data
    time = sf.getTimeBase('vrot')
    rho_pol = kk.fastRz2rho(31052,R,z,time)
    # compare with fast equilibrium
    diag = 'EQH'
    exp = 'EQXY'
    shot = 31052
    #Fix edition of the equilibrium
    out = kk.kkeqtp(shot,exp=exp,diag=diag)
    #delete timepoints outside of the equilibrium recontruction
    time = time
    rho = np.zeros([time.size,R.size])
    for jt in range(time.size):
        rho[jt] = kk.kkrzptfn(shot,time[jt],R,z,exp=exp,\
                            diag=diag,ed=out.ed).rho_p
    import dd
    dd = dd.shotfile('GQH',shot,experiment='AUGE')
    sig = dd('Raus')
    raus = sig.data
    raus_time = sig.time
    import matplotlib.pylab as plt
    plt.plot(time,(rho-rho_pol)*10.+2.1)
    plt.plot(raus_time,raus,'o')
    plt.show()
    ipsh()

def vesselStructure(nshot,ax=None):
    from matplotlib.patches import Polygon
    import matplotlib.pylab as plt
    show = False
    if ax == None:
        f = plt.figure()
        ax = f.add_subplot(111)
        show = True
    kk = KK()
    out = kk.kkgcread(nshot)
    for structure in out:
        if (structure['z'] < -1.3).any():
            #avoid to plot vessel
            continue
        ax.add_patch(Polygon(zip(structure['R']-0.001,structure['z']),facecolor='gray',edgecolor='none'))
    if show:
        ax.set_xlim([0,3])
        ax.set_ylim([-2,2])
        plt.show()


if __name__=="__main__":
    import matplotlib.pylab as plt
    kk = KK()
    a = kk.kkeqtpsa(36171,3.0,diag='EQI')
