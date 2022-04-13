from netCDF4 import Dataset
import parse_nml
import shutil, os, sys
import ctypes as ct
import numpy as np
from sf_dics import *

if ct.sizeof(ct.c_long) == 8:
  libso = 'lib64'
else:
  libso = 'lib'
libsfh = ct.cdll.LoadLibrary('/afs/ipp/aug/ads/'+libso+'/@sys/libsfh8.so')
#libsfh = ct.cdll.LoadLibrary('/afs/ipp/aug/ads/'+libso+'/@sys/libddww8.so.8.1.20111007')

class sfhhelp:
  status = False

class SFH:

  def Open(self,fname):
    """Opens the shot file header name, reads it in a temporary memory buffer and records the last modification date."""

    c_name = ct.c_char_p(fname)
    self.c_sfid = ct.c_int32(0)
    _sfid = ct.byref(self.c_sfid)

    err = libsfh.sfhopen(c_name,_sfid)
    libsfh.sfherror(err,'Open')
    return err

  def Close(self):
    """This routine does the following things:
- checks if sfid is still valid
- checks if the original shot file header is unchanged, otherwise there will be an error
- copies the original shot file header to name.BAK
- writes the new shot file header."""

    err=libsfh.sfhclose(self.c_sfid)
    libsfh.sfherror(err,'Close')

    return err

  def Lparrec(self,ps):
    """Returns a list of parameter record names, formats and number of values. Warning: listlen is used both as input and as output parameter. As input parameter it should contain the length of the lists. On return it contains the number of record names in this object."""

    c_ps=ct.c_char_p(ps)
    name_len=9
    listlen=100
    c_listlen=ct.c_uint32(listlen)
    _listlen=ct.byref(c_listlen)
    c_reclist = ((ct.c_char * name_len) * listlen)()
    c_rfmt=(ct.c_uint16 * listlen)()
    _rfmt = ct.byref(c_rfmt)
#    c_nrecs=(ct.c_uint32 * listlen)()
    c_nrecs=(ct.c_uint16 * listlen)()
    _nrecs = ct.byref(c_nrecs)

    err=libsfh.sfhlparrec(self.c_sfid,c_ps,_listlen,c_reclist,_rfmt,_nrecs)
    libsfh.sfherror(err,'Lparrec')

    output = sfhhelp()
    output.err = err

    nlist=c_listlen.value

    output.parrec={}

    for j in range(nlist):
      par_name = c_reclist[j][0:name_len].replace('\x00','').strip()
      fmt = c_rfmt[j]
      nrecs = c_nrecs[j]
      output.parrec[par_name]=(fmt,nrecs)
      print par_name,fmt,nrecs
    return output

  def Readparset(self,ps):
    """Reads all parameters of the object ps."""

    lpar  = self.Lparrec(ps)
    par_dic={}

    for pn in lpar.parrec.iterkeys():
      par_dic[pn]={}
      fmt   = lpar.parrec[pn][0]
      typ = fmt2type[fmt]
      pnlen = lpar.parrec[pn][1]
# Input
      c_ps=ct.c_char_p(ps)
      c_pn=ct.c_char_p(pn)
      c_typ=ct.c_uint32(typ)
      c_len=ct.c_uint32(pnlen)
# Output
      c_npar=ct.c_long(0)
      _npar=ct.byref(c_npar)

      if fmt in fmt2len.iterkeys():
        name_len=fmt2len[fmt]
        mylen = name_len * pnlen
        c_len=ct.c_uint32(mylen)
        c_data=(ct.c_char * mylen)()
        _data=c_data
        err=libsfh.sfhreadpar(self.c_sfid,c_ps,c_pn,c_typ,c_len,c_data,_npar)
        libsfh.sfherror(err,'readpar')
        out_arr=[''] * pnlen
        for j in range(pnlen):
          out_arr[j]=c_data[j*name_len:(j+1)*name_len].replace('\x00','')
      else:
        out_arr=np.empty(pnlen,dtype=type2np[typ])
        err=libsfh.sfhreadpar(self.c_sfid,c_ps,c_pn,c_typ,c_len, \
            out_arr.ctypes.data_as(ct.POINTER(ct.c_char)),_npar)
        libsfh.sfherror(err,'readpar')
      print pn,out_arr
      par_dic[pn]['data'] = out_arr
      par_dic[pn]['err'] = err
      par_dic[pn]['typ'] = typ

    return par_dic

  def Rdlist(self,obj):
    """Reads the members of a List to objnames. As a return value listlen shows the number of object names read."""

# Input
    c_obj=ct.c_char_p(obj)
# Output
    name_len=9
    listlen=900
    c_listlen=ct.c_long(listlen)
    _listlen=ct.byref(c_listlen)
    dlen = name_len * listlen
    clen=ct.c_long(dlen)
    c_objnames = ((ct.c_char * name_len) * listlen)()

    output = sfhhelp()

    err=libsfh.sfhrdlist(self.c_sfid,c_obj,_listlen,c_objnames)
    libsfh.sfherror(err,'rdlist')
    output.err = err

    if err != 0:
      nlist=c_listlen.value
      obj_name=[''] * nlist
      for j in range(nlist):
        obj_name[j] = c_objnames[j][0:name_len].replace('\x00','')
      output.objname = obj_name

    return output

  def Lonam(self):
    """Reads the object names and object types of a SFH to objnames and objtypes. Only listlen object names and object types will be read; the lists must be long enough to get all the names and types. As a return value listlen shows the number of objects in the SFH."""

# Output
    name_len=9
    listlen=600
    c_listlen=ct.c_uint32(listlen)
    _listlen=ct.byref(c_listlen)
    dlen = name_len * listlen
    clen=ct.c_long(dlen)
    _dlen=ct.byref(clen)
    c_objnames = ((ct.c_char * name_len) * listlen)()
    c_objtypes= (ct.c_uint16 * listlen)()
    _objtypes = ct.byref(c_objtypes)

    output = sfhhelp()

    err=libsfh.sfhlonam(self.c_sfid,_dlen,c_objnames,c_objtypes)
    libsfh.sfherror(err,'Lonam')
    output.err = err

    if err == 0:
      nlist=clen.value
      obj_name=['']
      output.lonam={}
      for j in range(nlist):
        obj_name = c_objnames[j][0:name_len].replace('\x00','')
        obj_type = c_objtypes[j]
        output.lonam[obj_name] = obj_type
    return output
    
  def Rdnsteps(self,obj):
    """Reads the number of time steps of a device to nsteps."""

# Input
    c_obj=ct.c_char_p(obj)
# Output
    c_nsteps=ct.c_long(0)
    _nsteps = ct.byref(c_nsteps)
    output = sfhhelp()

    err=libsfh.sfhrdnsteps(self.c_sfid,c_obj_nsteps)
    libsfh.sfherror(err,'rdnsteps')
    output.err = err

    output.nsteps={}
    if err == 0:
      output.nsteps[obj]=c_nsteps.value
    return output

  def Rdobj(self,obj):
    """Reads information of an object."""

# Input
    c_obj=ct.c_char_p(obj)
# Output
    c_objtyp=ct.uint16(0) #ct.c_long(0)
    _objtyp = ct.byref(c_objtyp)
    c_numdim=ct.c_long(0)
    _numdim = ct.byref(c_numdim)
    c_nsteps=ct.c_long(0)
    _nsteps = ct.byref(c_nsteps)
    c_format=ct.c_long(0)
    _format=ct.byref(c_format)
    output = sfhhelp()

    err=libsfh.sfhrdobj(self.c_sfid,c_obj,_objtyp,_numdim,_nsteps,_format)
    libsfh.sfherror(err,'rdobj')
    output.err = err
    output.rdobj={}
    if err == 0:
      output.rdobj['typ']   =c_objtyp.value
      output.rdobj['ndim']  =c_numdim.value
      output.rdobj['nsteps']=c_nsteps.value
      output.rdobj['format']=c_format.value
    return output

  def Rdindex24(self, signal_name):
    """Reads the higher indices of a signalgroup to index2, index3 and index4."""

# Input
    c_signame=ct.c_char_p(signal_name)
# Output
    c_index2=ct.c_uint32(0)
    _index2=ct.byref(c_index2)
    c_index3=ct.c_uint32(0)
    _index3=ct.byref(c_index3)
    c_index4=ct.c_uint32(0)
    _index4=ct.byref(c_index4)
    output = sfhhelp()
    
    err=libsfh.sfhrdindex24(self.c_sfid,c_signame,_index2,_index3,_index4)
    libsfh.sfherror(err,'rdmap')
    output.err = err
    output.rdindex24={}
    if err == 0:
      output.rdindex24[2] = c_index2.value
      output.rdindex24[3] = c_index3.value
      output.rdindex24[4] = c_index4.value
    return output

  def Rdarea(self,obj):
    """Reads the number of time steps and the three dimensions of an areabase."""

# Input
    c_obj=ct.c_char_p(obj)
# Output
    c_nsteps=ct.c_long(0)
    _nsteps = ct.byref(c_nsteps)
    c_x=ct.c_long(0)
    _x=ct.byref(c_x)
    c_y=ct.c_long(0)
    _y=ct.byref(c_y)
    c_z=ct.c_long(0)
    _z=ct.byref(c_z)
    output = sfhhelp()

    err=libsfh.sfhrdarea(self.c_sfid,c_obj,_nsteps,_x,_y,_z)
    libsfh.sfherror(err,'rdarea')
    output.err = err
    output.rdarea={}
    if err == 0:
      output.rdarea['nsteps']=c_nsteps.value
      output.rdarea['nx']=c_x.value
      output.rdarea['ny']=c_y.value
      output.rdarea['nz']=c_z.value
 
    return output

  def Rdmap(self,signal_name,index2,index3,index4):
    """Search all Mapping tables for an appropriate record."""

# Input
    c_signame=ct.c_char_p(signal_name)
    c_index2=ct.c_uint32(index2)
    c_index3=ct.c_uint32(index3)
    c_index4=ct.c_uint32(index4)
# Output
    name_len=9
    c_devname=(ct.c_char * name_len)()
    c_chan=ct.c_uint16(0)
    _chan=ct.byref(c_chan)
    output = sfhhelp()

    err=libsfh.sfhrdmap(self.c_sfid,c_signame,c_index2,c_index3,c_index4,c_devname, _chan)
    libsfh.sfherror(err,'rdmap')
    output.err = err
    output.rdmap={}
    if err == 0:
      output.rdmap['devname'] = c_devname.value
      output.rdmap['chan'] = c_chan.value

    return output

  def Relations(self,obj):
    """Creates a list rels with the names of all relations defined and the number of relations in numrels. """

# FIXME: does not work

    print 'Begin ',obj
# Input
    c_obj=ct.c_char_p(obj)
# Output
    name_len=9
    nrel=2
    c_rels = ((ct.c_char * name_len) * nrel)()
    _rels=ct.byref(c_rels)
    c_nrel=ct.c_uint32(nrel)
    _nrel=ct.byref(c_nrel)
    output = sfhhelp()

    err=libsfh.sfhrelations(self.c_sfid,c_obj,c_rels,_nrel)
    libsfh.sfherror(err,'Relations')
    output.err = err

    if err == 0:
      nlist=c_nrel.value
      print nlist
      obj_name=[''] * nlist
      for j in range(nlist):
        obj_name[j] = c_rels[j][0:name_len].replace('\x00','')
      output.relations = obj_name

    return output

#===============
#  Modification
#===============

  def Mdarea(self,obj,nsteps,nx,ny,nz):

    c_obj=ct.c_char_p(obj)
    c_nsteps=ct.c_uint32(nsteps)
    c_nx=ct.c_uint32(nx)
    c_ny=ct.c_uint32(ny)
    c_nz=ct.c_uint32(nz)

    err=libsfh.sfhmdarea(self.c_sfid,c_obj,c_nsteps,c_nx,c_ny,c_nz)
    libsfh.sfherror(err,'mdarea')

  def Mdindex(self,obj,nx,ny,nz,qual=False):

    c_obj=ct.c_char_p(obj)
    c_nx=ct.c_uint32(nx)
    c_ny=ct.c_uint32(ny)
    c_nz=ct.c_uint32(nz)

    if qual:
      err=libsfh.sfhmdqualindex(self.c_sfid,c_obj,c_nx,c_ny,c_nz)
    else:
      err=libsfh.sfhmdindex24(self.c_sfid,c_obj,c_nx,c_ny,c_nz)
    libsfh.sfherror(err,'mdindex24')

  def Modtim(self,obj,nt):

    c_obj=ct.c_char_p(obj)
    c_nt=ct.c_uint32(nt)

    err=libsfh.sfhmodtim(self.c_sfid,c_obj,c_nt)
    libsfh.sfherror(err,'modtim')

  def Modindex1(self,obj,nt):

    c_obj=ct.c_char_p(obj)
    c_nt=ct.c_uint32(nt)

    err=libsfh.sfhmdindex1(self.c_sfid,c_obj,c_nt)
    libsfh.sfherror(err,'modindex1')

  def Modsgr(self,obj,dims,qual=False):

    nx=1+np.zeros(4,dtype=int)
    for jdim in range(len(dims)):
      nx[jdim]=dims[jdim]
    self.Modtim(obj,nx[0])
    self.Mdindex(obj,nx[1],nx[2],nx[3],qual=qual)

  def Modpar(self,ps,pn,dat):

    lpar  = self.Lparrec(ps)
    par_dic={}

    par_dic[pn]={}
    fmt   = lpar.parrec[pn][0]
    typ = fmt2type[fmt]
    pnlen = len(dat)
# Input
    c_ps=ct.c_char_p(ps)
    c_pn=ct.c_char_p(pn)
    c_typ=ct.c_uint32(typ)
    c_len=ct.c_uint32(pnlen)
    _typ=ct.byref(c_typ)
    _len=ct.byref(c_len)
# Output

    print 'Modpar test:',pn,typ,pnlen,len(dat)
    if fmt in fmt2len.iterkeys():
      name_len=fmt2len[fmt]
      c_len = ct.c_uint32(pnlen)
      c_data=ct.c_char_p(dat)
      err=libsfh.sfhmodpar(self.c_sfid,c_ps,c_pn,c_typ,c_len,c_data)
    else:
      c_data = (ct.c_long * pnlen)()
      for j in range(pnlen):
        c_data[j]=ct.c_long(dat[j])
      _data=ct.byref(c_data)

      err=libsfh.sfhmodpar(self.c_sfid,c_ps,c_pn,c_typ,c_len,dat.ctypes.data_as(ct.POINTER(ct.c_char)))
    libsfh.sfherror(err,'Modpar')

#===============
#  New Methods
#===============

  def Mapping(self):
    
    objects = self.Lonam().lonam
    mapping = []

    for i in objects:
      if obj_d[objects[i]] == "Sig_Group":
        indices = self.Rdindex24(i).rdindex24
        for i2 in range(indices[2]):
          for i3 in range(indices[3]):
            for i4 in range(indices[4]):
              mapping.append([i, [i2, i3, i4], 
                            self.Rdmap(i, i2 + 1 ,i3 + 1, i4 + 1).rdmap['chan']])

    return mapping

  def Channelmapping(self):
  
    raw_mapping = self.Mapping()
    mapping_table = {}
    for i in range(len(raw_mapping)):
      mapping_table[raw_mapping[i][2]] = raw_mapping[i][0] + '_' + str(raw_mapping[i][1][0]) + '-' + str(raw_mapping[i][1][1]) + '-' + str(raw_mapping[i][1][2])
    
    mapping_array=[]
    for i in range(len(raw_mapping)):
      mapping_array.append(mapping_table[i])
    
    return np.array(mapping_array)

  def Signalmapping(self):
  
    raw_mapping = self.Mapping()
    mapping_table = {}
    for i in range(len(raw_mapping)):
      mapping_table[raw_mapping[i][0] + '_' + str(raw_mapping[i][1][0]) + '-' + str(raw_mapping[i][1][1]) + '-' + str(raw_mapping[i][1][2])] = raw_mapping[i][2]

    return mapping_table

  def sfhmod(self, cdf_file, nml='', fsfh='TRA00000.sfh', source='/afs/ipp/home/t/transp/pub/TRA00000.sfh.temp'):

    if not os.path.isfile(cdf_file):
      print(cdf_file+' not found')
      sys.exit()
    runid = cdf_file.split('/')[-1][0:8]
    cdf   = Dataset(cdf_file, 'r', format='NETCDF4')
    sfh_d={}
    mydic={}
    
    try:
      shutil.copy2(source, fsfh)
    except:
      print('Could not copy file '+source+' to '+fsfh)
    err= self.Open(fsfh)

    if err != 0:
      return
    lnm = self.Lonam()
    if lnm.err != 0:
      return

    for obj,devtyp in lnm.lonam.iteritems():
      if devtyp == 8:
        sfh_d[obj]={}
        sfh_d[obj]['devtyp']=devtyp
        sfh_d[obj]['data']=cdf.variables[obj][:]
        sfh_d[obj]['dim']=cdf.variables[obj].dimensions
        nt_cdf=sfh_d[obj]['data'].shape[0]
        self.Modtim(obj,nt_cdf)

    for obj,devtyp in lnm.lonam.iteritems():

# Devtyp: 13 Areabase, 6 SignalGroup, 7 Signal, 8 TimeBase, 4 ParameterSet

# Parameter Sets

      if devtyp == 4:
        if nml == '': # Read values from namelist
          print('Parameter sets not read from any namelist')
        else:
          sfh_d[obj] = {}
          sfh_d[obj]['devtyp'] = devtyp
          pars = self.Lparrec(obj)
          sfh_d[obj]['info'] = pars.parrec
          for pn, val in pars.parrec.iteritems():
            sfh_d[obj][pn]={}
            if pn == 'runid':
              sfh_d[obj][pn]['data'] = runid
              sfh_d[obj][pn]['format'] = 1794
            else:
              parse = parse_nml.PARSE_NML(nml, pn, fmt=val[0])
              sfh_d[obj][pn]['data'] = parse.outarr
              sfh_d[obj][pn]['format'] = val[0]
            print pn, sfh_d[obj][pn]['data']

# AB, TB, SIG, SGR

      else:
        name_flg=False
        if obj in cdf.variables.iterkeys():
          name_flg=True
          cdfobj=obj
        else:
          if len(obj) == 8:
            name_flg=False
            for cdfobj in cdf.variables.iterkeys():
              if obj in cdfobj:
                print 'CDF signal '+cdfobj+' stored as shotfile '+obj
                name_flg=True
                break
        if not name_flg:
          print 'Signal '+obj+' not found in '+cdf_file
        else:
          sfh_d[obj]={}
          sfh_d[obj]['devtyp']=devtyp
          sfh_d[obj]['data']=cdf.variables[cdfobj][:]
          sfh_d[obj]['dim']=cdf.variables[cdfobj].dimensions
          nt_cdf=sfh_d[obj]['data'].shape[0]
          if len(sfh_d[obj]['data'].shape) ==  2:
            nx_cdf=sfh_d[obj]['data'].shape[1]
          if devtyp == 13:
            if len(sfh_d[obj]['data'].shape) ==  1:
              nt_cdf=1
              nx_cdf=sfh_d[obj]['data'].shape[0]
            self.Mdarea(obj,1,nx_cdf,0,0) # wwainsert not yet ready for AB(t)
#            self.Mdarea(obj,nt_cdf,nx_cdf,0,0)
          if devtyp == 6:
            self.Mdindex(obj,nx_cdf,0,0)
          if devtyp in (6,7):
            self.Modtim(obj,nt_cdf)


    self.Close()
    return sfh_d
