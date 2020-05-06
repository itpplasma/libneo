""" 
Write/Read data in UFILES format (version 0.2)

UFILES is an ASCII format defined in the manual:
  http://w3.pppl.gov/~pshare/help/ufiles_manual.html

How to use ufiles.py:
  
  import ufiles as uf

======
Read:
======

  u = uf.RU(fname)

  Input:
    fname [string] is the filename with the full path 

  Output:
    u.shot: shot number as appear in the first line of the ufiles
    u.dim: rank of the function, i.e. rnk(f(X)) = 1, rnk(F(X,Y))=2 
    u.label['X','Y','Z']: dict of the independent variables' labels
    u.values['X','Y','Z']: dict of numpy arrays for the X, Y and Z grid
    u.flabel: label of the dependent variable
    u.fvalues: numpy array of the function data

  Example:
    plt.plot(u.values['Y'],u.fvalues[0,:])

======
Write:
======

  u = uf.WU(u_dic,udir=path)

  Input:
    path [string] is the target dir; default is $HOME/udb/$shot
    u_dic is the input dictionary:
      u_dic['pre']  = output file prefix [string]
      u_dic['ext']  = output file extension [string]
      u_dic['shot'] = shot number [string]
      u_dic['scal'] [list]
        u_dic['scal']=[[label1 [string],scal_val1 [float]], [label2, scal_val2], ...]
      u_dic['X'] [dictionary]
        u_dic['X']['lbl'] = grid variable label + unit [string]
        u_dic['X']['arr'] = grid variable array [numpy 1D array]
      Same for 'Y' or 'Z' in case of 2D, 3D ufiles
      u_dic['data'] [dictionary]
        u_dic['data']['lbl'] = independent variable label + unit [string]
        u_dic['data']['arr'] = independent variable array [numpy 1D, 2D or 3D array]


"""

__author__  = 'Giovanni Tardini (Tel. 1898)'
__version__ = '0.02'
__date__    = '13.03.2012'

import datetime,os,sys,shutil
import numpy as np
import rw_for

now = datetime.datetime.now()
w_for = rw_for.RW_FOR().wr_for
sspl  = rw_for.RW_FOR().ssplit

class WU:

    def __init__(self, uf_d, udir=''):

        sshot  = uf_d['shot']
        prestr = uf_d['pre']
        extstr = uf_d['ext']
        comment = 'Generated with write_u.py'
        if 'comm' in uf_d.iterkeys():
            comment = uf_d['comm']
        if udir == '':
            udir = '%s/udb/%s' %(os.getenv('HOME'), sshot)
        os.system('mkdir -p %s' %udir)
        uf = '%s/%s%s.%s' %(udir, prestr, sshot, extstr)
        if os.path.isfile(uf):
            ufbak = '%s~' %uf
            print('cp %s %s' %(uf, ufbak))
            shutil.copy2(uf, ufbak)
        print('Writing file %s' %uf)

        dev = 'AUGD'
        xyz = ('X', 'Y', 'Z')
        if 'grid' in uf_d.iterkeys():
            ndim = len(uf_d['grid'])
        else:
            ndim=0

# Header

        ufs = '  %5s%4s %1i 0 6              ;-SHOT #- F(X) DATA -UF%1iDWR- %s\n'\
             %(sshot,  dev,  ndim,  ndim, now.strftime('%d%b%Y'))
        ufs += ''.ljust(30) + ';-SHOT DATE-  UFILES ASCII FILE SYSTEM\n'

        if 'scal' in uf_d.iterkeys():
            nscal = len(uf_d['scal'])
            ufs += (' %2d' %nscal).ljust(30) + \
                ';-NUMBER OF ASSOCIATED SCALAR QUANTITIES-\n'
            for jscal in range(nscal):
                lbl, val = uf_d['scal'][jscal]
                ufs += ('%11.4E' %val).ljust(30) + ';-SCALAR,  LABEL FOLLOWS:\n'
                ufs += lbl.ljust(30) + '\n'
        else:
            ufs += '  0'.ljust(30) +';-NUMBER OF ASSOCIATED SCALAR QUANTITIES-\n'

        if ndim > 0:
            if (ndim == 1):
                if 'unit' in uf_d['grid']['X'].iterkeys():
                    ufs += '%-20s%-10s;-INDEPENDENT VARIABLE LABEL-\n' \
                    %(uf_d['grid']['X']['name'], uf_d['grid']['X']['unit'])
                else:
                    ufs += '%-30s;-INDEPENDENT VARIABLE LABEL-\n' %uf_d['grid']['X']['lbl']
            else:
                for coord in xyz:
                    if coord in uf_d['grid'].iterkeys():
                        if 'unit' in uf_d['grid'][coord].iterkeys():
                            ufs += '%-20s%-10s;-INDEPENDENT VARIABLE LABEL: ' \
                            %(uf_d['grid'][coord]['name'], uf_d['grid'][coord]['unit'])
                        else:
                            ufs += '%-30s;-INDEPENDENT VARIABLE LABEL: ' %uf_d['grid'][coord]['lbl']
                        ufs += '%1s-\n' %coord

            if 'unit' in uf_d['data'].iterkeys():
                ufs += '%-20s%-10s;-DEPENDENT VARIABLE LABEL-\n' \
                %(uf_d['data']['name'], uf_d['data']['unit'])
            else:
                ufs += '%-30s;-DEPENDENT VARIABLE LABEL-\n'   %uf_d['data']['lbl']
            ufs += '0'.ljust(30) + ';-PROC CODE- 0:RAW 1:AVG 2:SM. 3:AVG+SM\n'

# Grids

            for coord in xyz:
                if coord in uf_d['grid'].iterkeys():
                    nl = len(uf_d['grid'][coord]['arr'])
                    ufs += ('   %7d' %nl).ljust(30) + ';-# OF %1s PTS-\n' %coord
            for coord in xyz:
                if coord in uf_d['grid'].iterkeys():
                    ufs += w_for(uf_d['grid'][coord]['arr'])

# Data

            ufs += w_for(uf_d['data']['arr'])

        ufs += ';----END-OF-DATA-----------------COMMENTS:-----------\n'
        ufs = ufs.replace('\n', '\n ')
        ufs += comment + '\n'

        f = open(uf, 'w')
        f.write(ufs)
        f.close()


class RU:


    def __init__(self, fname):

        nldata = 6 # N. of values per line

# === Open the file
        f = open(fname, 'r')
        self.file = fname
        self.lines = f.readlines()
        self.label = {}
        self.values = {}
        self.n = {}
        self.vars = ['X', 'Y', 'Z']

# === Headers
        self.nvar = 0
        nlines = len(self.lines)

        for j in range(nlines):
            l = self.lines[j]
            a = l.split()

            if l.find('-SHOT #') !=-1 or l.find('; Shot #') != -1:
                self.shot = a[0][:5]
                self.dim = int(a[1])

            if l.find('-INDEPENDENT VARIABLE LABEL') !=-1 or l.find('; Independent variable') != -1:
                for var in self.vars:
                    if l.find(var + '-') != -1:
                        self.label[var] = l.split(';-')[0][1: ]
                        self.nvar +=1
                if self.nvar == 0:
                    self.label['X'] = l.split(';-')[0][1:] if '-INDEPENDENT' in l else l.split('; Independent')[0][1:]
                    self.nvar = 1

            for var in self.vars:
                if l.find('-# OF %s PTS-' %var[0]) != -1:
                    self.n[var] = int(a[0])
                if ';-# of radial pts  %s'%var in l:
                    self.n[var] = int(a[0])

            if l.find('-DEPENDENT VARIABLE LABEL-') != -1: 
                self.flabel = l.split(';-')[0][1: ]

# === Data
        if self.dim != self.nvar:
            print('Inconsistency in the number of independent variables dim=%i nvar=%i'%(self.dim, self.nvar))
            sys.exit()

        list_var = self.vars[ :self.nvar]

        jstart = 3 + 2*self.dim + 2
        jlast = jstart

        while True: # skip over lines until we find first true data line
            try:
                temp = np.array(self.lines[jlast].split(), dtype=float)
                del temp
                break
            except:
                jlast += 1

        for var in list_var:
            self.values[var] = np.zeros(self.n[var])
            jstart = jlast
            jlast = self.n[var]/nldata + jstart
            if self.n[var] %nldata != 0:
                jlast += 1
            i1 = 0

            for j in range(jstart, jlast):
                l = self.lines[j]
                a = l.split()
                for l in range(len(a)):
                    self.values[var][i1+l] = float(a[l])
                i1 += len(a)
            self.nf = self.n[var]

# ... F
        if self.dim == 1:
            ftmp = np.zeros(self.n['X'])
        elif self.dim == 2:
            ftmp = np.zeros(self.n['X']*self.n['Y'])  
        elif self.dim == 3:
            ftmp = np.zeros(self.n['X']*self.n['Y']*self.n['Z'])

        jstart = jlast
        jlast = len(ftmp)/nldata + jstart
        if len(ftmp) %nldata != 0:
                jlast += 1 
        i1 = 0
        jn = 0
        for j in range(jstart, jlast):
            l = self.lines[j]
            a = l.split()
            if l.find('END-OF-DATA') !=-1:
                break
            a = sspl(l)
            i2 = i1 + len(a)
            ftmp[i1: i2] = a
            i1 = i1 + nldata
            jn += nldata

        if self.dim == 1:
            self.fvalues = ftmp

        elif self.dim == 2:
            self.fvalues = ftmp.reshape(self.n['Y'], self.n['X']).T

        elif self.dim == 3:
            self.fvalues = ftmp.reshape(self.n['Z'], self.n['Y'], self.n['X']).T
