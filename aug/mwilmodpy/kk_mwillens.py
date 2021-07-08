import ctypes as ct
import numpy as np
from IPython import embed
import os
# Input: byref for arrays and strings, ctypes for scalars
# Output: always byref


if ct.sizeof(ct.c_long) == 8:
    libso = 'lib64'
else:
    libso = 'lib'
kklib = '/afs/ipp/aug/ads/'+libso+'/@sys/libkk8.so'

if not os.path.isfile(kklib):
    kklib = '/afs/ipp/aug/ads/'+libso+'/amd64_sles15/libkk8.so'
libkk = ct.cdll.LoadLibrary(kklib)


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

    def kkeqrzcl(self, nshot,tshot, exp=exp_eq, diag=dia_eq, ed=0):
        
        nt = np.size(tshot)

        # Input
        c_ed    = ct.c_long(ed)
        _ed     = ct.byref(c_ed)
        c_err   = ct.c_long(0)
        _err    = ct.byref(c_err)
        c_shot  = ct.c_long(nshot)
        _shot   = ct.byref(c_shot)
        c_exp   = ct.c_char_p(exp)
        c_dia   = ct.c_char_p(diag)

        
 #       ncl=13
        # Output
#        c_ncl  = ct.c_long(ncl)
#        _ncl   = ct.byref(c_ncl)
        pncl=40
        c_ncl  = ct.c_int(pncl)
        _ncl   = ct.byref(c_ncl)

        rcl   = (ct.c_float * pncl)()
        zcl   = (ct.c_float * pncl)()
        cld  = (ct.c_float * pncl)()
        cle  = (ct.c_float * pncl)()
        _rcl  = ct.byref(rcl)
        _zcl  = ct.byref(zcl)
        _cld = ct.byref(cld)
        _cle = ct.byref(cle)



# Numpy ououtput     = kkhelp()tput

        output     = kkhelp()
        output.ncl = np.zeros((1))
        output.rcl = np.zeros((nt))
        output.zcl = np.zeros((nt))
        output.cld = np.zeros((nt))
        output.cle = np.zeros((nt))
        embed()
        for jt, t in enumerate(np.ravel(tshot)):
            c_t    = ct.c_float(t)
            _tshot = ct.byref(c_t)
            
            status = libkk.kkeqrzcl(_err, c_exp, c_dia, c_shot, _ed, _tshot,
                                    _ncl, _rcl,_zcl,_cld, _cle)

            
            output.ncl = c_ncl.value
            output.rcl[jt] = rcl
            output.zcl[jt] = zcl
            output.cld[jt] = cld
            output.cle[jt] = cle
        
            
        output.ed  = c_ed.value
        output.err = c_err.value
#    Numpy output
        nt = c_ntim.value
        output.time = np.array(tim_kk[0:nt])
        return output
        



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
        Rdim = np.size(np.shape(Rin))
        if Rdim == 1:
            nr = len(Rin)
        else:
            nr = len(Rin[0])
       
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
        if Rdim == 1:
            c_Rin    = (ct.c_float * nr)()
            c_Rin[:] = Rin[:]
            _Rin     = ct.byref(c_Rin)
            c_zin    = (ct.c_float * nr)()
            c_zin[:] = zin[:]
            _zin     = ct.byref(c_zin)
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
            c_t    = ct.c_float(t)
            _tshot = ct.byref(c_t)
            if Rdim == 2:
                c_Rin    = (ct.c_float * nr)()
                c_Rin[:] = Rin[jt,:]
                _Rin     = ct.byref(c_Rin)
                c_zin    = (ct.c_float * nr)()
                c_zin[:] = zin[jt,:]
                _zin     = ct.byref(c_zin)

            status = libkk.kkrzptfn(_err, c_exp, c_dia, c_shot, _ed, _tshot, \
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
            c_t    = ct.c_float(t)
            _tshot = ct.byref(c_t)

            status = libkk.kkrhopto(_err, c_exp, c_dia, c_shot, _ed, _tshot, \
                                    _rhop, c_nr, _rhot, _fpf, _ftf)
            output.rho_t[jt, :] = rhot[:]
            output.pf   [jt, :] = fpf[:]
            output.tf   [jt, :] = ftf[:]

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



    def kkrhorz(self, nshot, tshot_in, rhop, angle=0, exp=exp_eq, diag=dia_eq, ed=0):

        nt=np.size(tshot_in)
        nDimRhop = np.size(np.shape(rhop))
        if ((nDimRhop==2) & (np.shape(rhop)[0]==nt)):
            dimRhop = True
            nr=np.size(rhop[0])
        else:
            nr = np.size(rhop)
            dimRhop = False
            
            #nr = np.size(rhop)
   #     nr = np.size(tshot)
        nDimAngle = np.size(np.shape(angle))
        nt=np.size(tshot_in)
        
        if nDimAngle == 1:
            nAngle = np.size(angle)
            angle = np.ravel(angle)
        elif nDimAngle == 2:
            idxNr = np.where(np.array(np.shape(angle))==nr)[0]
            if np.size(idxNr) == 0:
                angle = np.ravel(angle)
                nAngle = np.size(angle)
            elif np.size(idxNr) == 1:
                angle = np.swapaxes(angle,idxNr,1)
                nAngle = np.size(angle[:,0])
            elif np.size(idxNr) == 2:
                angle = angle
                nAngle = nr
            else:
                angle = np.ravel(angle)
                nAngle = nr
        else:
            angle = np.ravel(angle)
            nAngle = 1            

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


# Special input
        c_iorg    = ct.c_long(0)
        _iorg     = ct.byref(c_iorg)
        c_rhop    = (ct.c_float * nr)()
        if dimRhop==False:
            c_rhop[:] = rhop[:]
            _rhop     = ct.byref(c_rhop)
            
        lexp      = ct.c_long( len(exp) )
        ldia      = ct.c_long( len(diag) )

# Output
        c_zn = (ct.c_float * nr)()
        c_rn = (ct.c_float * nr)()
        _zn = ct.byref(c_zn)
        _rn = ct.byref(c_rn)

# Numpy output
        output = kkhelp()
        
        output.r = np.zeros((nt,nAngle, nr))
        output.z    = np.zeros((nt,nAngle, nr))

        
        if nDimAngle != 2:
            for t,tshot in enumerate(tshot_in): 
                for i,an in enumerate(np.ravel(angle)):
                    c_angle   = ct.c_float(an)
                    _angle    = ct.byref(c_angle)
                    c_tshot = ct.c_float(tshot)
                    _tshot  = ct.byref(c_tshot)
                    if dimRhop==True:
                        c_rhop[:] = rhop[t,:]
                        _rhop     = ct.byref(c_rhop)
                        
                    status = libkk.kkrhorz(_err, c_exp, c_dia, c_shot, _ed, _tshot,  \
                               _rhop, c_nr, c_iorg, c_angle, _rn, _zn)
 

                    output.r[t,i, :] = np.array(c_rn[0:nr])
                    output.z[t,i, :] = np.array(c_zn[0:nr])

        else:

            c_nr    = ct.c_long(1)
            _nr     = ct.byref(c_nr)

            c_rhop    = (ct.c_float * 1)() 
            _rhop     = ct.byref(c_rhop)

            c_zn = (ct.c_float * 1)()
            c_rn = (ct.c_float * 1)()
            _zn = ct.byref(c_zn)
            _rn = ct.byref(c_rn)
            
            for t,tshot in enumerate(tshot_in): 
                for j in (np.arange(nr)):
                    for i in (np.arange(nAngle)):

                        c_tshot = ct.c_float(tshot)
                        _tshot  = ct.byref(c_tshot)
                    
                        c_angle   = ct.c_float(angle[i,j])
                        if dimRhop==True:
                            c_rhop[:] = [rhop[t,j]]

                        status = libkk.kkrhorz(_err, c_exp, c_dia, c_shot, _ed, _tshot,  \
                               _rhop, c_nr, c_iorg, c_angle, _rn, _zn)
                                        
                        output.r[t,i, j] = np.array(c_rn[0])
                        output.z[t,i, j] = np.array(c_zn[0])

                        output.ed  = c_ed.value
                        output.err = c_err.value

# Numpy output
        if nAngle == 1:
            output.r = np.squeeze(output.r)
            output.z = np.squeeze(output.z) 

        return output



#kkrhoToP (iERR  ,expnam,dianam,nSHOT,nEDIT,tSHOT,
#C               >                 rhoTF,Lrho,
#C               <                 rhoPF,      fPF, fTF)

    def kkrhotop(self, nshot, tshot, rhot, exp=exp_eq, diag=dia_eq, ed=0):

        nr = np.size(rhot)
   #     nr = np.size(tshot)
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
        c_rhot    = (ct.c_float * nr)()
        c_rhot[:] = rhot[:]
        _rhot     = ct.byref(c_rhot)
        lexp      = ct.c_long( len(exp) )
        ldia      = ct.c_long( len(diag) )

# Output
        c_rhop = (ct.c_float * nr)()     
        _rhop = ct.byref(c_rhop)
        c_fpf = (ct.c_float * nr)()
        _fpf = ct.byref(c_fpf)
        c_ftf = (ct.c_float * nr)()
        _ftf = ct.byref(c_ftf)


        status = libkk.kkrhotop(_err, c_exp, c_dia, c_shot, _ed, _tshot,  \
                               _rhot, c_nr, _rhop, _fpf, _ftf)
 

        output = kkhelp()
        output.rhop  = np.array(c_rhop[0:nr])
        output.fpf  = np.array(c_fpf[0:nr])
        output.ftf  = np.array(c_ftf[0:nr])
        output.ed  = c_ed.value
        output.err = c_err.value
        

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

    def kkrzBrzt(self, nshot, tshot_in, Rin, zin, exp=exp_eq, diag=dia_eq, ed=0):

# Input
        nr = len(Rin)
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


# Special input
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

        
        output = kkhelp()
        
        nt=np.size(tshot_in)
        if nt>1:
            
            output = kkhelp()
            output.br  = np.zeros((nt, nr)) 
            output.bz  = np.zeros((nt, nr)) 
            output.bt  = np.zeros((nt, nr)) 
            output.fpf = np.zeros((nt, nr)) 
            output.fjp = np.zeros((nt, nr))
            for jt, tshot in enumerate(np.ravel(tshot_in)):
                c_tshot = ct.c_float(tshot)
                _tshot  = ct.byref(c_tshot)
                status = libkk.kkrzbrzt(_err, c_exp, c_dia, c_shot, _ed, _tshot, _Rin, _zin, \
                                c_nr, _br, _bz, _bt, _fpf, _fjp)
                
                output.ed  = c_ed.value
                output.err = c_err.value

                # Numpy output
                output.br[jt]  = np.array(br[0:nr])
                output.bz[jt]  = np.array(bz[0:nr])
                output.bt[jt]  = np.array(bt[0:nr])
                output.fpf[jt] = np.array(fpf[0:nr])
                output.fjp[jt] = np.array(fjp[0:nr])
    
        else:
            c_tshot = ct.c_float(tshot_in)
            _tshot  = ct.byref(c_tshot)

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



    def kkrztmprzt(self, nshot, tshot, dt, Rin, zin, pin, exp=exp_eq, diag='YGC', ed=0):

# Input
        nr = len(Rin)
        c_nr    = ct.c_long(nr)
        _nr     = ct.byref(c_nr)
        c_ed    = ct.c_long(ed)
        _ed     = ct.byref(c_ed)
        c_err   = ct.c_long(0)
        _err    = ct.byref(c_err)
        c_shot  = ct.c_long(nshot)
        _shot   = ct.byref(c_shot)
        c_exp   = ct.c_char_p(exp)
        _exp   = ct.byref(c_exp)
        c_dia   = ct.c_char_p(diag)
        _dia   = ct.byref(c_dia)
        c_tshot = ct.c_float(tshot)
        _tshot  = ct.byref(c_tshot)
        c_dt = ct.c_float(dt)
        _dt  = ct.byref(c_dt)

# Special input
        c_Rin = (ct.c_float * nr)()
        c_Rin[:] = Rin[:]
        _Rin = ct.byref(c_Rin)
        c_zin = (ct.c_float * nr)()
        c_zin[:] = zin[:]
        _zin = ct.byref(c_zin)
        c_pin = (ct.c_float * nr)()
        c_pin[:] = pin[:]
        _pin = ct.byref(c_pin)

# Output mpr, mpz, mpp
        mpr = (ct.c_float * nr)()
        mpz = (ct.c_float * nr)()
        mpp  = (ct.c_float * nr)()
        _mpr  = ct.byref(mpr)
        _mpz  = ct.byref(mpz)
        _mpp  = ct.byref(mpp)
        
        print 'tshot; ',tshot

        status = libkk.kkrztmprzt(_err, _exp, _dia, c_shot, _ed,c_tshot, c_dt, \
                                      c_Rin, c_zin, c_pin, c_nr, \
                                      _mpr, _mpz, _mpp )

#libkk.kkrztmprzt(byref(error), byref(exp), byref(diag), shot, byref(edit), time, dt,
#                 rwant, zwant, pwant, lwant,
#                 byref(mpr), byref(mpz), byref(mpp) ) 

        output = kkhelp()
        output.ed  = c_ed.value
        output.err = c_err.value

# Numpy output
        output.mpr  = np.array(mpr[0:nr])
        output.mpz  = np.array(mpz[0:nr])
        output.mpt  = np.array(mpp[0:nr])


        return output


    def kkrzq(self, nshot, tshot, Rin, zin, exp=exp_eq, diag=dia_eq, ed=0):

# Input
        Rdim = np.size(np.shape(Rin))
        if Rdim == 1:
            nr = len(Rin)
        else:
            nr = len(Rin[0])
       
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
        if Rdim == 1:
            c_Rin    = (ct.c_float * nr)()
            c_Rin[:] = Rin[:]
            _Rin     = ct.byref(c_Rin)
            c_zin    = (ct.c_float * nr)()
            c_zin[:] = zin[:]
            _zin     = ct.byref(c_zin)
        elif Rdim == 2:
            c_Rin    = (ct.c_float * nr)()
            c_Rin[:] = Rin[0]
            _Rin     = ct.byref(c_Rin)
            c_zin    = (ct.c_float * nr)()
            c_zin[:] = zin[0]
            _zin     = ct.byref(c_zin)
         
# Output
        fpf   = (ct.c_float * nr)()
        q   = (ct.c_float * nr)()
        rhop  = (ct.c_float * nr)()
        
        _fpf  = ct.byref(fpf)
        _q  = ct.byref(q)
        _rhop = ct.byref(rhop)
        

        output = kkhelp()

# Numpy output
        output.pf    = np.zeros((nt, nr))
        output.q    = np.zeros((nt, nr))
        output.rho_p = np.zeros((nt, nr))


        for jt, t in enumerate(np.ravel(tshot)):
            c_t    = ct.c_float(t)
            _tshot = ct.byref(c_t)
            if Rdim == 2:
                c_Rin    = (ct.c_float * nr)()
                c_Rin[:] = Rin[jt,:]
                _Rin     = ct.byref(c_Rin)
                c_zin    = (ct.c_float * nr)()
                c_zin[:] = zin[jt,:]
                _zin     = ct.byref(c_zin)
                status = libkk.kkrzq(_err, c_exp, c_dia, c_shot, _ed, _tshot, \
                                      c_Rin, c_zin, c_nr, \
                                      _q, _fpf, _rhop )
            elif Rdim == 1:
                status = libkk.kkrzq(_err, c_exp, c_dia, c_shot, _ed, _tshot, \
                                      c_Rin, c_zin, c_nr, \
                                      _q, _fpf, _rhop )     
                
                
            output.pf[jt, :] = fpf[:]
            output.q[jt, :] = q[:]
            output.rho_p[jt, :] = rhop[:]


        output.ed    = c_ed.value
        output.err   = c_err.value
        output.pf    = np.squeeze(output.pf)
        output.q    = np.squeeze(output.q)
        output.rho_p = np.squeeze(output.rho_p)
        

        return output

#C     ___________________ calculate from rhoPF-values  -->   q   _values
#C                kkrhoPFq (iERR  ,expnam,dianam,nSHOT,nEDIT,tSHOT,
#C               >                 rhoPF,Lrho,
#C               <                 fq   ,      fPF)


    def kkrhopfq(self, nshot, tshot, rhop, exp=exp_eq, diag=dia_eq, ed=0):

        nr = np.size(rhop)
   #     nr = np.size(tshot)
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
        c_rhop    = (ct.c_float * nr)()
        c_rhop[:] = rhop[:]
        _rhop     = ct.byref(c_rhop)
        lexp      = ct.c_long( len(exp) )
        ldia      = ct.c_long( len(diag) )

# Output
        c_q = (ct.c_float * nr)()     
        _q = ct.byref(c_q)
        c_fpf = (ct.c_float * nr)()
        _fpf = ct.byref(c_fpf)
 
        status = libkk.kkrhopfq(_err, c_exp, c_dia, c_shot, _ed, _tshot,  \
                               _rhop, c_nr, _q, _fpf )
 
        output = kkhelp()
        output.q  = np.array(c_q[0:nr])
        output.fpf  = np.array(c_fpf[0:nr])
        output.ed  = c_ed.value
        output.err = c_err.value
        

        return output





    def kkrhotpq(self, nshot, tshot, rhot, exp=exp_eq, diag=dia_eq, ed=0):

        nr = np.size(rhot)
   #     nr = np.size(tshot)
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
        c_rhot    = (ct.c_float * nr)()
        c_rhot[:] = rhot[:]
        _rhot     = ct.byref(c_rhot)
        lexp      = ct.c_long( len(exp) )
        ldia      = ct.c_long( len(diag) )

# Output
        c_q = (ct.c_float * nr)()     
        _q = ct.byref(c_q)
        c_fpf = (ct.c_float * nr)()
        _fpf = ct.byref(c_fpf)
        c_ftf = (ct.c_float * nr)()
        _ftf = ct.byref(c_ftf)

 
        status = libkk.kkrhotpq(_err, c_exp, c_dia, c_shot, _ed, _tshot,  \
                               _rhot, c_nr, _q, _fpf, _ftf )
 
        output = kkhelp()
        output.q  = np.array(c_q[0:nr])
        output.fpf  = np.array(c_fpf[0:nr])
        output.ftf  = np.array(c_ftf[0:nr])
        output.ed  = c_ed.value
        output.err = c_err.value
        

        return output



