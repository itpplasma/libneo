import ctypes as ct
import numpy as np

# Input: byref for arrays and strings, ctypes for scalars
# Output: always byref

if ct.sizeof(ct.c_long) == 8:
    libso = 'lib64'
else:
    libso = 'lib'
libkk = ct.cdll.LoadLibrary('/afs/ipp/aug/ads/'+libso+'/@sys/libkk.so')


exp_eq='AUGD'
dia_eq='EQH'

class kkhelp:
    status = False

class KK:

    def kkeqints(self,isw):

# Input
        c_err = ct.c_long(0)
        _err  = ct.byref(c_err)
        c_isw = ct.c_long(isw)
        _isw  = ct.byref(c_isw)

        status = libkk.kkeqints(_err,c_isw)
        return status

    def kkeqtp(self, nshot, exp=exp_eq, diag=dia_eq, ed=0):

        ntim = 20000
# Input
        c_ed    = ct.c_long(ed)
        _ed     = ct.byref(c_ed)
        c_err   = ct.c_long(0)
        _err    = ct.byref(c_err)
        c_shot  = ct.c_long(nshot)
        _shot   = ct.byref(c_shot)
        c_exp   = ct.c_char_p(exp)
        c_dia   = ct.c_char_p(diag)
# Special input
        c_nprnt = ct.c_long(0)
        _nprnt  = ct.byref(c_nprnt)
# Output
        c_ntim  = ct.c_int(ntim)
        _ntim   = ct.byref(c_ntim)
        tim_kk  = (ct.c_float * ntim)()
        _tim_kk = ct.byref(tim_kk)

        status = libkk.kkeqtp(_err, c_exp, c_dia, c_shot, _ed, _nprnt, _ntim, _tim_kk)

        output     = kkhelp()
        output.ed  = c_ed.value
        output.err = c_err.value
#    Numpy output
        nt = c_ntim.value
        output.time = np.array(tim_kk[0:nt])
        return output

    def kkrzptfn(self, nshot, tshot, Rin, zin, exp=exp_eq, diag=dia_eq, ed=0):

# Input
        nr = len(Rin)
        nt = np.size(tshot)
        c_nr   = ct.c_long(nr)
        _nr    = ct.byref(c_nr)
        c_ed   = ct.c_long(ed)
        _ed    = ct.byref(c_ed)
        c_err  = ct.c_long(0)
        _err   = ct.byref(c_err)
        c_shot = ct.c_long(nshot)
        _shot  = ct.byref(c_shot)
        c_exp  = ct.c_char_p(exp)
        c_dia  = ct.c_char_p(diag)
        
# Special input
        c_Rin = (ct.c_float * nr)()
        c_Rin[:] = Rin[:]
        _Rin = ct.byref(c_Rin)
        c_zin = (ct.c_float * nr)()
        c_zin[:] = zin[:]
        _zin = ct.byref(c_zin)
# Output
        fpf   = (ct.c_float * nr)()
        ftf   = (ct.c_float * nr)()
        rhop  = (ct.c_float * nr)()
        rhot  = (ct.c_float * nr)()
        _fpf  = ct.byref(fpf)
        _ftf  = ct.byref(ftf)
        _rhop = ct.byref(rhop)
        _rhot = ct.byref(rhot)

        output = kkhelp()

# Numpy output
        output.pf    = np.zeros((nt, nr))
        output.tf    = np.zeros((nt, nr))
        output.rho_p = np.zeros((nt, nr))
        output.rho_t = np.zeros((nt, nr))

        for jt, t in enumerate(np.ravel(tshot)):
            c_t = ct.c_float(t)
            _time = ct.byref(c_t)
            status = libkk.kkrzptfn(_err, c_exp, c_dia, c_shot, _ed, _time, \
                                    _Rin, _zin, c_nr, _fpf, _rhop, _ftf, _rhot)
                
            output.pf   [jt, :] = fpf[:]
            output.tf   [jt, :] = ftf[:]
            output.rho_p[jt, :] = rhop[:]
            output.rho_t[jt, :] = rhot[:]

        output.ed    = c_ed.value
        output.err   = c_err.value
        output.pf    = np.squeeze(output.pf)
        output.tf    = np.squeeze(output.tf)
        output.rho_p = np.squeeze(output.rho_p)
        output.rho_t = np.squeeze(output.rho_t)

        return output

    def kkrhopto(self, nshot, tshot, rhop, exp=exp_eq, diag=dia_eq, ed=0):

        nr = len(rhop)
        nt = np.size(tshot)

# Input
        c_nr   = ct.c_long(nr)
        _nr    = ct.byref(c_nr)
        c_ed   = ct.c_long(ed)
        _ed    = ct.byref(c_ed)
        c_err  = ct.c_long(0)
        _err   = ct.byref(c_err)
        c_shot = ct.c_long(nshot)
        _shot  = ct.byref(c_shot)
        c_exp  = ct.c_char_p(exp)
        c_dia  = ct.c_char_p(diag)
# Special input
        c_rhop    = (ct.c_float * nr)()
        c_rhop[:] = rhop[:]
        _rhop     = ct.byref(c_rhop)
        lexp      = ct.c_long( len(exp) )
        ldia      = ct.c_long( len(diag) )
# Output
        fpf   = (ct.c_float * nr)()
        ftf   = (ct.c_float * nr)()
        rhot  = (ct.c_float * nr)()
        _fpf  = ct.byref(fpf)
        _ftf  = ct.byref(ftf)
        _rhot = ct.byref(rhot)
        
        output = kkhelp()
 
# Numpy output
        output.rho_t = np.zeros((nt, nr))
        output.pf    = np.zeros((nt, nr))
        output.tf    = np.zeros((nt, nr))

        for jt, t in enumerate(np.ravel(tshot)):
            c_t   = ct.c_float(t)
            _time = ct.byref(c_t)

            status = libkk.kkrhopto(_err, c_exp, c_dia, c_shot, _ed, _time, \
                                    _rhop, c_nr, _rhot, _fpf, _ftf)
            output.rho_t[jt, :] = rhot[:]
            output.pf[jt, :]    = fpf[:]
            output.tf[jt, :]    = ftf[:]

        output.ed    = c_ed.value
        output.err   = c_err.value
        output.pf    = np.squeeze(output.pf)
        output.tf    = np.squeeze(output.tf)
        output.rho_t = np.squeeze(output.rho_t)

        return output

    def kkeqpsp(self, nshot, tshot, flux, exp=exp_eq, diag=dia_eq, ed=0):

# Input
        c_ed    = ct.c_long(ed)
        c_err   = ct.c_long(0)
        _err    = ct.byref(c_err)
        _ed     = ct.byref(c_ed)
        c_shot  = ct.c_long(nshot)
        _shot   = ct.byref(c_shot)
        c_exp   = ct.c_char_p(exp)
        c_dia   = ct.c_char_p(diag)
        c_tshot = ct.c_float(tshot)
        _tshot  = ct.byref(c_tshot)
        c_tsfh  = ct.c_float(0.)
        _tsfh   = ct.byref(c_tsfh)
        lexp    = ct.c_long( len(exp) )
        ldia    = ct.c_long( len(diag) )
# Specific input
        c_flux = ct.c_float(flux)
        _flux  = ct.byref(c_flux)
# Output
        nrzmax = 1000
        npfm   = 1001
        mpfm   = 1001
        pfm    = (ct.c_float * mpfm * npfm )()
        _pfm   = ct.byref(pfm)
        c_ndim = ct.c_long(1000)
        _ndim  = ct.byref(c_ndim)
        rz     = (ct.c_float * nrzmax * 2)()
        _rz    = ct.byref(rz)
        c_nrz  = ct.c_long(nrzmax)
        _nrz   = ct.byref(c_nrz)
        status = libkk.kkeqpsp(_err, c_exp, c_dia, c_shot, _ed, _tshot,  \
                               c_flux, _ndim, _rz, _nrz, _pfm, _tsfh)

        print ''
        print 'time = ', c_tsfh.value
        output     = kkhelp()
        output.ed  = c_ed.value
        output.err = c_err.value
# Numpy output
        n_rz = c_nrz.value
        if (n_rz > 0):
            output.r_surf = np.array(rz[0][0:n_rz])
            output.z_surf = np.array(rz[1][0:n_rz])
        else:
            output.r_surf = np.zeros(1)
            output.z_surf = np.zeros(1)
        return output

    def kkeqtfl(self, nshot, tshot, exp=exp_eq, diag=dia_eq, ed=0):

        ndxymax = 3200
# Input
        c_ed    = ct.c_long(ed)
        _ed     = ct.byref(c_ed)
        c_err   = ct.c_long(0)
        _err    = ct.byref(c_err)
        c_shot  = ct.c_long(nshot)
        _shot   = ct.byref(c_shot)
        c_exp   = ct.c_char_p(exp)
        c_dia   = ct.c_char_p(diag)
        c_tshot = ct.c_float(tshot)
        _tshot  = ct.byref(c_tshot)
# Output
        c_pfl = (ct.c_float * ndxymax)()
        c_tfl = (ct.c_float * ndxymax)()
        c_lpf = ct.c_long(ndxymax)
        _pfl  = ct.byref(c_pfl)
        _tfl  = ct.byref(c_tfl)
        _lpf  = ct.byref(c_lpf)

        status = libkk.kkeqtfl(_err, c_exp, c_dia, c_shot, _ed, _tshot,  \
                               _lpf, _pfl, _tfl)

        output     = kkhelp()
        output.ed  = c_ed.value
        output.err = c_err.value
# Numpy output
        if c_err.value != 0:
            print ' ERROR: kkeqtfl', exp, diag
        else:
            nr = c_lpf.value+1
            output.pf = np.array(c_pfl[0:nr])
            output.tf = np.array(c_tfl[0:nr])
        return output

    def kkpfrhop(self, nshot, tshot, fpf, exp=exp_eq, diag=dia_eq, ed=0):

        nr = len(fpf)
# Input
        c_ed    = ct.c_long(ed)
        _ed     = ct.byref(c_ed)
        c_err   = ct.c_long(0)
        _err    = ct.byref(c_err)
        c_shot  = ct.c_long(nshot)
        _shot   = ct.byref(c_shot)
        c_exp   = ct.c_char_p(exp)
        c_dia   = ct.c_char_p(diag)
        c_tshot = ct.c_float(tshot)
        _tshot  = ct.byref(c_tshot)
# Special input
        c_fp    = (ct.c_float * nr)()
        c_fp[:] = fpf[:]
        _fp     = ct.byref(c_fp)
        c_nr    = ct.c_long(nr)
        _nr     = ct.byref(c_nr)
# Output
        c_rhop = (ct.c_float * nr)()
        _rhop  = ct.byref(c_rhop)

        status = libkk.kkpfrhop(_err, c_exp, c_dia, c_shot, _ed, _tshot,  \
                                _fp, c_nr, _rhop)

        output     = kkhelp()
        output.ed  = c_ed.value
        output.err = c_err.value
# Numpy output
        output.rho_p = np.array(c_rhop[0:nr])
        return output

    def kktfrhot(self, nshot, tshot, ftf, exp=exp_eq, diag=dia_eq, ed=0):

        nr = len(ftf)
# Input
        c_ed    = ct.c_long(ed)
        _ed     = ct.byref(c_ed)
        c_err   = ct.c_long(0)
        _err    = ct.byref(c_err)
        c_shot  = ct.c_long(nshot)
        _shot   = ct.byref(c_shot)
        c_exp   = ct.c_char_p(exp)
        c_dia   = ct.c_char_p(diag)
        c_tshot = ct.c_float(tshot)
        _tshot  = ct.byref(c_tshot)
# Special input
        c_ft    =    (ct.c_float * nr)()
        c_ft[:] = ftf[:]
        _ft     = ct.byref(c_ft)
        c_nr    = ct.c_long(nr)
        _nr     = ct.byref(c_nr)
# Output
        c_rhot = (ct.c_float * nr)()
        _rhot  = ct.byref(c_rhot)

        status = libkk.kktfrhot(_err, c_exp, c_dia, c_shot, _ed, _tshot,  \
                                _ft, c_nr, _rhot)
        
        output     = kkhelp()
        output.ed  = c_ed.value
        output.err = c_err.value
# Numpy output
        output.rho_t = np.array(c_rhot[0:nr])
        return output

    def kkrhotpq(self, nshot, tshot, rhot, exp=exp_eq, diag=dia_eq, ed=0):

        nr = len(rhot)
# Input
        c_ed    = ct.c_long(ed)
        _ed     = ct.byref(c_ed)
        c_err   = ct.c_long(0)
        _err    = ct.byref(c_err)
        c_shot  = ct.c_long(nshot)
        _shot   = ct.byref(c_shot)
        c_exp   = ct.c_char_p(exp)
        c_dia   = ct.c_char_p(diag)
        c_tshot = ct.c_float(tshot)
        _tshot  = ct.byref(c_tshot)
# Special input
        c_nr      = ct.c_long(nr)
        _nr       = ct.byref(c_nr)
        c_rhot    = (ct.c_float * nr)()
        c_rhot[:] = rhot[:]
        _rhot     = ct.byref(c_rhot)
        lexp      = ct.c_long( len(exp) )
        ldia      = ct.c_long( len(diag) )
# Output
        c_qp = (ct.c_float * nr)()
        c_fp = (ct.c_float * nr)()
        c_ft = (ct.c_float * nr)()
        _qp  = ct.byref(c_qp)
        _fp  = ct.byref(c_fp)
        _ft  = ct.byref(c_ft)

        status = libkk.kkrhotpq(_err, c_exp, c_dia, c_shot, _ed, _tshot,  \
                                _rhot, c_nr, _qp, _fp, _ft)

        output     = kkhelp()
        output.ed  = c_ed.value
        output.err = c_err.value
# Numpy output
        output.fp = np.array(c_fp[0:nr])
        output.ft = np.array(c_ft[0:nr])
        output.qp = np.array(c_qp[0:nr])
        return output

    def kkeqpfx(self, nshot, tshot, exp = exp_eq, diag = dia_eq, ed = 0):

        nr = 4
        nr1 = nr+1
# Input
        c_ed    = ct.c_long(ed)
        _ed     = ct.byref(c_ed)
        c_err   = ct.c_long(0)
        _err    = ct.byref(c_err)
        c_shot  = ct.c_long(nshot)
        _shot   = ct.byref(c_shot)
        c_exp   = ct.c_char_p(exp)
        c_dia   = ct.c_char_p(diag)
        c_tshot = ct.c_float(tshot)
        _tshot  = ct.byref(c_tshot)
# Special input
        c_lpfx = ct.c_long(nr)
        _lpfx  = ct.byref(c_lpfx)
# Output
        c_pfx  = (ct.c_float * nr1)()
        c_rpfx = (ct.c_float * nr1)()
        c_zpfx = (ct.c_float * nr1)()
        _pfx   = ct.byref(c_pfx)
        _rpfx  = ct.byref(c_rpfx)
        _zpfx  = ct.byref(c_zpfx)

        status = libkk.kkeqpfx(_err, c_exp, c_dia, c_shot, _ed, _tshot,  \
                               _lpfx, _pfx, _rpfx, _zpfx)

        output     = kkhelp()
        output.ed  = c_ed.value
        output.err = c_err.value
# Float output
        output.Raxis = c_rpfx[0]
        output.zaxis = c_zpfx[0]
        output.Rspx  = c_rpfx[1]
        output.zspx  = c_zpfx[1]
        output.Rlim  = c_rpfx[2]
        output.zlim  = c_zpfx[2]
        return output

    def kkeqpres(self, nshot, tshot, exp=exp_eq, diag=dia_eq, ed=0):

        print 'Reading pressure profile'
        npres_max = 100000
# Input
        c_ed    = ct.c_long(ed)
        _ed     = ct.byref(c_ed)
        c_err   = ct.c_long(0)
        _err    = ct.byref(c_err)
        c_shot  = ct.c_long(nshot)
        _shot   = ct.byref(c_shot)
        c_exp   = ct.c_char_p(exp)
        c_dia   = ct.c_char_p(diag)
        c_tshot = ct.c_float(tshot)
        _tshot  = ct.byref(c_tshot)
# Output
        c_npres = ct.c_long(npres_max)
        _npres  = ct.byref(c_npres)
        c_fp    = (ct.c_float * npres_max)()
        c_pres  = (ct.c_float * npres_max)()
        c_presp = (ct.c_float * npres_max)()
        _fp     = ct.byref(c_fp)
        _pres   = ct.byref(c_pres)
        _presp  = ct.byref(c_presp)

        status = libkk.kkeqpres(_err, c_exp, c_dia, c_shot, _ed, _tshot,  \
                                _npres, _fp, _pres, _presp)

        n_pres     = c_npres.value+1
        output     = kkhelp()
        output.ed  = c_ed.value
        output.err = c_err.value
# Numpy output
        output.pres = np.array(c_pres[0:n_pres])
        output.fp   = np.array(c_fp[0:n_pres])
        return output

    def kkeqqpl(self, nshot, tshot, exp=exp_eq, diag=dia_eq, ed=0):

        lpf_max = 100000
# Input
        c_ed    = ct.c_long(ed)
        _ed     = ct.byref(c_ed)
        c_err   = ct.c_long(0)
        _err    = ct.byref(c_err)
        c_shot  = ct.c_long(nshot)
        _shot   = ct.byref(c_shot)
        c_exp   = ct.c_char_p(exp)
        c_dia   = ct.c_char_p(diag)
        c_tshot = ct.c_float(tshot)
        _tshot  = ct.byref(c_tshot)
# Output
        c_lpf = ct.c_long(lpf_max)
        _lpf  = ct.byref(c_lpf)
        c_pfl = (ct.c_float * lpf_max)()
        c_qpl = (ct.c_float * lpf_max)()
        _pfl  = ct.byref(c_pfl)
        _qpl  = ct.byref(c_qpl)

        status = libkk.kkeqqpl(_err, c_exp, c_dia, c_shot, _ed, _tshot,  \
                               _lpf, _pfl, _qpl)

        output = kkhelp()
        output.ed  = c_ed.value
        output.err = c_err.value
#    Numpy output
        nrho = c_lpf.value+1
        output.pfl = np.array(c_pfl[0:nrho])
        output.qpl = np.array(c_qpl[0:nrho])
        return output

    def kkeqffs(self, nshot, tshot, exp=exp_eq, diag=dia_eq, ed=0):

        nr_max = 100000
# Input
        c_ed    = ct.c_long(ed)
        _ed     = ct.byref(c_ed)
        c_err   = ct.c_long(0)
        _err    = ct.byref(c_err)
        c_shot  = ct.c_long(nshot)
        _shot   = ct.byref(c_shot)
        c_exp   = ct.c_char_p(exp)
        c_dia   = ct.c_char_p(diag)
        c_tshot = ct.c_float(tshot)
        _tshot  = ct.byref(c_tshot)
# Special input
        c_typ = ct.c_long(11)
        _typ  = ct.byref(c_typ)
# Output
        c_nr     = ct.c_long(nr_max)
        _nr      = ct.byref(c_nr)
        c_fp     = (ct.c_float * nr_max)()
        c_rbphi  = (ct.c_float * nr_max)()
        c_drbphi = (ct.c_float * nr_max)()
        _fp      = ct.byref(c_fp)
        _rbphi   = ct.byref(c_rbphi)
        _drbphi  = ct.byref(c_drbphi)

        status = libkk.kkeqffs(_err, c_exp, c_dia, c_shot, _ed, _tshot,  \
                               c_typ, _nr, _fp, _rbphi, _drbphi)

        n_rbphi = c_nr.value+1
        output = kkhelp()
        output.ed  = c_ed.value
        output.err = c_err.value
# Numpy output
        output.fp = np.array(c_fp[0:n_rbphi])
        output.rb = np.array(c_rbphi[0:n_rbphi])
        return output

    def kkcutsx(self, nshot, tshot, Rin, zin, phi, exp=exp_eq, diag=dia_eq, ed=0):

# Input
        c_ed    = ct.c_long(ed)
        _ed     = ct.byref(c_ed)
        c_err   = ct.c_long(0)
        _err    = ct.byref(c_err)
        c_shot  = ct.c_long(nshot)
        _shot   = ct.byref(c_shot)
        c_exp   = ct.c_char_p(exp)
        c_dia   = ct.c_char_p(diag)
        c_tshot = ct.c_float(tshot)
        _tshot  = ct.byref(c_tshot)
# Special input
        c_Rin = ct.c_float(Rin)
        c_zin = ct.c_float(zin)
        c_phi = ct.c_float(phi)
        _Rin  = ct.byref(c_Rin)
        _zin  = ct.byref(c_zin)
        _phi  = ct.byref(c_phi)
#Output
        c_Rout  = ct.c_float(0)
        c_zout  = ct.c_float(0)
        c_psisx = ct.c_float(0)
        c_tsfh  = ct.c_float(0)
        _Rout   = ct.byref(c_Rout)
        _zout   = ct.byref(c_zout)
        _psisx  = ct.byref(c_psisx)
        _tsfh   = ct.byref(c_tsfh)

        status = libkk.kkcutsx(_err, c_exp, c_dia, c_shot, _ed, _tshot,  \
                               c_Rin, c_zin, c_phi, _Rout, _zout, _psisx, _tsfh)

        output = kkhelp()
        output.ed  = c_ed.value
        output.err = c_err.value
        output.Rsx = c_Rout.value
        output.zsx = c_zout.value
        return output

    def kkrhorz(self, nshot, tshot, rhop, angle=0, exp=exp_eq, diag=dia_eq, ed=0):

        nr = np.size(rhop) 
        na = np.size(angle) 

# Input
        c_nr    = ct.c_long(nr)
        _nr     = ct.byref(c_nr)
        c_ed    = ct.c_long(ed)
        _ed     = ct.byref(c_ed)
        c_err   = ct.c_long(0)
        _err    = ct.byref(c_err)
        c_shot  = ct.c_long(nshot)
        _shot   = ct.byref(c_shot)
        c_exp   = ct.c_char_p(exp)
        c_dia   = ct.c_char_p(diag)
        c_tshot = ct.c_float(tshot)
        _tshot  = ct.byref(c_tshot)
# Special input
        c_iorg    = ct.c_long(0)
        _iorg     = ct.byref(c_iorg)
        c_rhop    = (ct.c_float * nr)()
        c_rhop[:] = rhop
        _rhop     = ct.byref(c_rhop)
        lexp      = ct.c_long( len(exp) )
        ldia      = ct.c_long( len(diag) )
        
        
        angle = np.atleast_1d(angle)
        R = np.empty((na,nr),dtype='float32')
        Z = np.empty((na,nr),dtype='float32')
        
        for ia, a in enumerate(angle):
            c_angle   = ct.c_float(a)
            _angle    = ct.byref(c_angle)
    # Output
            c_zn = (ct.c_float * nr)()
            c_rn = (ct.c_float * nr)()
            _zn = ct.byref(c_zn)
            _rn = ct.byref(c_rn)

            status = libkk.kkrhorz(_err, c_exp, c_dia, c_shot, _ed, _tshot,  \
                                _rhop, c_nr, c_iorg, c_angle, _rn, _zn)
            R[ia] = c_rn[:nr]
            Z[ia] = c_zn[:nr]

            

        output = kkhelp()
        output.ed  = c_ed.value
        output.err = c_err.value
# Numpy output
        output.r = np.squeeze(R)  #squalze to make it backward compatible
        output.z = np.squeeze(Z)

        return output

    def kkGCd0(self, nshot = 28746, ndim = 800, ndgc = 45, exp = exp_eq, diag = 'YGC', ed = 0):

# Input
        c_ed   = ct.c_long(ed)
        c_err  = ct.c_long(0)
        c_shot = ct.c_long(nshot)
        c_exp  = ct.c_char_p(exp)
        c_dia  = ct.c_char_p(diag)
        _ed    = ct.byref(c_ed)
        _err   = ct.byref(c_err)
        _shot  = ct.byref(c_shot)
# Special input
        c_ndim = ct.c_long(ndim)
        c_ndgc = ct.c_long(ndgc)
        _ndim  = ct.byref(c_ndim)
        _ndgc  = ct.byref(c_ndgc)
# Output
        c_xygc  = ((ct.c_float * ndim) * 2)()
        c_ngc   = ct.c_long(0)
        c_ixbeg = (ct.c_int * ndgc)()
        c_lenix = (ct.c_int * ndgc)()
        c_valix = (ct.c_int * ndgc)()
        c_gcnam = (ct.c_char_p * ndgc)()
        _xygc   = ct.byref(c_xygc)
        _ngc    = ct.byref(c_ngc)
        _ixbeg  = ct.byref(c_ixbeg)
        _lenix  = ct.byref(c_lenix)
        _valix  = ct.byref(c_valix)
        _gcnam  = ct.byref(c_gcnam)

        status = libkk.kkgcd0(_err, c_exp, c_dia, c_shot, _ed, c_ndim,  \
                              _xygc, c_ndgc, _ngc, _ixbeg, _lenix, _valix, _gcnam)
        
        output = kkhelp()
        output.ed  = c_ed.value
        output.err = c_err.value

        output.gc_x = {}
        output.gc_y = {}
        print ndgc
        for jgc in range(ndgc):
            if c_valix[jgc] > 0:
                jl   = c_ixbeg[jgc] - 1
                jlen = c_lenix[jgc]
                lbl  = c_gcnam[jgc][0:8]
                output.gc_x[lbl] = np.array(c_xygc[0][jl:jl+jlen])
                output.gc_y[lbl] = np.array(c_xygc[1][jl:jl+jlen])
                
        return output

    def kkrzBrzt(self, nshot, tshot, Rin, zin, exp=exp_eq, diag=dia_eq, ed=0):

# Input
        nr = np.size(Rin)


        c_nr    = ct.c_long(nr)
        _nr     = ct.byref(c_nr)
        c_ed    = ct.c_long(ed)
        _ed     = ct.byref(c_ed)
        c_err   = ct.c_long(0)
        _err    = ct.byref(c_err)
        c_shot  = ct.c_long(nshot)
        _shot   = ct.byref(c_shot)
        c_exp   = ct.c_char_p(exp)
        c_dia   = ct.c_char_p(diag)
        c_tshot = ct.c_float(tshot)
        _tshot  = ct.byref(c_tshot)

# Special input
        Rin = np.atleast_1d(Rin)
        zin = np.atleast_1d(zin)
        c_Rin = (ct.c_float * nr)()
        c_Rin[:] = Rin[:]
        _Rin = ct.byref(c_Rin)
        c_zin = (ct.c_float * nr)()
        c_zin[:] = zin[:]
        _zin = ct.byref(c_zin)

# Output    Br, Bz, Bt,      fPF, fJp
        fpf = (ct.c_float * nr)()
        fjp = (ct.c_float * nr)()
        br  = (ct.c_float * nr)()
        bz  = (ct.c_float * nr)()
        bt  = (ct.c_float * nr)()
        _fpf = ct.byref(fpf)
        _fjp = ct.byref(fjp)
        _br  = ct.byref(br)
        _bz  = ct.byref(bz)
        _bt  = ct.byref(bt)

        status = libkk.kkrzbrzt(_err, c_exp, c_dia, c_shot, _ed, _tshot, _Rin, _zin, \
                                c_nr, _br, _bz, _bt, _fpf, _fjp)

        output = kkhelp()
        output.ed  = c_ed.value
        output.err = c_err.value

# Numpy output
        output.br  = np.array(br[0:nr])
        output.bz  = np.array(bz[0:nr])
        output.bt  = np.array(bt[0:nr])
        output.fpf = np.array(fpf[0:nr])
        output.fjp = np.array(fjp[0:nr])

        return output
