#""" Module providing the means to read ASDEX Upgrade data
#
#Example:
#
#sf = dd.shotfile()
#if sf.Open(diagnostic,shotnumber):
#   timebase = sf.GetTimebase(signal)
#   signal = sf.GetSignal(signal)
#   sf.Close() 
#"""
import ctypes as ct
import numpy as np
#print "my dd"
# rdm 13.08.12: lbuf is ct.c_uint32, no longer ct.c_int32

# Source: /usr/ads/codes/formats
fmt_d = { \
  1:'BYTE', 2:'CHAR', 3:'SHORT_INT', 4:'INTEGER', 5:'IEEE_FLOAT', \
  6:'IEEE_DOUBLE', 7:'LOGICAL', 8:'CHAR_REAL', 9:'U_SHORT', 10:'IBM_REAL', \
  11:'IBM_DOUBLE', 12:'CONDENSED_TB', 13:'LONGLONG', 14:'U_INTEGER', \
  15:'U_LONGLONG', 1794:'CHAR_8', 3842:'CHAR16', 7938:'CHAR_32', \
  12034:'CHAR_48', 16130:'CHAR_64', 18178:'CHAR_72' \
}

# Source: http://www.aug.ipp.mpg.de/aug/manuals/ddww/libddww8.html#ddobjval
type_d = { \
  0:'RAW', 1:'INTEGER', 2:'REAL', 3:'DOUBLE PRECISION', 4: 'COMPLEX', \
  5:'LOGICAL', 6:'CHARACTER', 10:'LONG LONG', 11:'SHORT', 12:'UNSIGNED SHORT', \
  13:'UNSIGNED INTEGER', 14:'UNSIGNED LONG LONG' \
}

fmt2type = { 1:5, 2:6, 3:11, 4:1, 5:2, 6:3, 9:12, 14:13, 13:2,15:14, \
             1794:6, 3842:6, 7938:6, 12034:6, 16130:6, 18178:6}
fmt2len = {2:1, 1794:8, 3842:16, 7938:32, 12034:48, 16130:64, 18178:72}

type2ct = { 1:ct.c_int32, 2:ct.c_float, 3:ct.c_double, \
  6:ct.c_char, 11:ct.c_int16, 12:ct.c_uint16 }

fmt2ct = { 2:ct.c_char, 3:ct.c_int16, 4:ct.c_int32, 5:ct.c_float, \
   6:ct.c_double, 9:ct.c_uint16 }

obj_d = { 2:'LIST', 4:'PARSET', 5:'MAP', 6:'SGR', 7:'SIG', 8:'TB', \
          9:'SFLIST', 13:'AB', 14:'QUAL' }

#libddww = ct.cdll.LoadLibrary('/afs/ipp-garching.mpg.de/aug/ads/lib64/@sys/libddww8.so.8.1.20111007')
libddww = ct.cdll.LoadLibrary('/afs/ipp-garching.mpg.de/aug/ads/lib64/@sys/libddww8.so.8.1')

class dd_info:
   status = False

class shotfile:
   """ Shotfile class to read ASDEX Upgrade data """

   def __init__(self):
      self.__diaref = ct.c_int(0)
      self.__status = False
      self.units = None

   def __del__(self):
      self.Close()

   def Open(self, diagname, shotnumber, experiment='AUGD', edition=0):
      """
        Open a shotfile \n\n
      """

      if self.__status:
         self.Close()
      error     = ct.c_int(0)
      _error    = ct.byref(error)
      self.edit = ct.c_int(edition)
      _edition  = ct.byref(self.edit)
      cshot     = ct.c_uint32(0)
      _cshot    = ct.byref(cshot)
      shot      = ct.c_uint32(int(shotnumber))
      _shot     = ct.byref(shot)
      diag      = ct.c_char_p(diagname)
      exp       = ct.c_char_p(experiment)
      dat       = 18*'d'
      date      = ct.c_char_p(dat)   
      _diaref   = ct.byref(self.__diaref)
      lexp  = ct.c_ulonglong(len(experiment))
      ldiag = ct.c_ulonglong(len(diagname))
      ldate = ct.c_ulonglong(len(dat))

      #result = libddww.ddcshotnr_(exp, diag, _shot, _cshot, lexp, ldiag)

      #if result == 0 and cshot.value == int(shotnumber):
      result = libddww.ddopen_(_error, exp, diag, _shot, _edition, _diaref, \
                                  date, lexp, ldiag, ldate)
      if result == 0:
         self.__status = True
         self.edition = self.edit.value
      else:
         self.__diaref = ct.c_int(0)
         self.GetError(error)
      return self.__status

   def Close(self):
      """ Close a shotfile """
      if self.__status:
         error = ct.c_int(0)
         _error = ct.byref(error)
         _diaref = ct.byref(self.__diaref)

         result = libddww.ddclose_(_error,_diaref)

         self.GetError(error)   
         self.__status = False
         self.__diaref = ct.c_int(0)
         print 'DDclose: edition '+str(self.edition)+' closed'

   def GetError(self, error):
      """ Print error code """ 
      if error.value != 0:
         _error = ct.byref(error)
         ID    = ""
         _ID   = ct.c_char_p(ID)
         lID   = ct.c_ulonglong(len(ID))
         ctrl  = ct.c_uint32(3)
         _ctrl = ct.byref(ctrl)

         result = libddww.xxerror_(_error, _ctrl, _ID, lID)

         print ID

   def GetInfo(self, name):
      """ Returns information about the specified signal."""
# Assuming: not more than one TB, not more than one AB 
      if self.__status:

         output = dd_info
         rel = self.GetRelations(name)
         output.error = rel.error
         if rel.error == 0:
            output.error = 0
            output.tname = None
            output.aname = None
            jtime = None
            jarea = None
            output.tlen = None
            for jid, id in enumerate(rel.typ):
               if id == 8:
                  jtime = jid
                  output.tname = rel.txt[jid]
               if id == 13:
                  jarea = jid
                  output.aname = rel.txt[jid]
                  error   = ct.c_int32(0)
                  _error  = ct.byref(error)
                  _diaref = ct.byref(self.__diaref)
                  _name   = ct.c_char_p(name)
                  sizes   = (ct.c_uint32*3)()
                  _sizes  = ct.byref(sizes)
                  adim    = (ct.c_int32*3)()
                  _adim   = ct.byref(adim)
                  index   = ct.c_int32(0)
                  _index  = ct.byref(index)
                  lname   = ct.c_ulonglong(len(name))

                  result = libddww.ddainfo_(_error, _diaref, _name, _sizes, \
                                            _adim, _index, lname)
                  if error.value == 0:
                     output.index = index.value
                  else:
                     output.index = None

            head = self.GetObjectHeader(name)
            output.error = head.error
            if head.error == 0:
               output.buf    = head.buffer
               output.objtyp = output.buf[0]
               output.fmt    = output.buf[14]
               len_dim       = output.buf[16]
               dims          = np.array(output.buf[18:22][::-1], dtype=np.int32)
               output.ind    = np.array(dims[dims > 0])

               if output.objtyp == 8: # If 'name' is a TB
                  output.tlen = output.buf[21] # = dims[0]
                  output.tfmt = output.buf[14]
               else:
                  if jtime != None:
# Check consistency with TB length
                     thead = self.GetObjectHeader(rel.txt[jtime])
                     tbuf = thead.buffer
                     output.tlen = tbuf[21]
                     output.tfmt = tbuf[14]
                  else:
                     print('No TB found for %s %s' %(obj_d[output.objtyp], name))

               if output.objtyp == 13: # If 'name' is an AB
                  output.atlen = output.buf[21]
                  output.afmt  = output.buf[14]
                  sizes = np.array(output.buf[18:21], dtype = np.int32)
                  output.sizes = sizes[sizes > 0]
               else:
# Beware: other than in DDAINFO2, here 'sizes' can have less than 
# 3 dims, as the 0-sized are removed. Usually (always?) it has 1 dim.
                  if jarea != None:
                     ahead = self.GetObjectHeader(rel.txt[jarea])
                     abuf = ahead.buffer
                     output.atlen = abuf[21] # #time points of AB
                     output.afmt  = abuf[14]
                     sizes = np.array(abuf[18:21], dtype = np.int32)
                     output.sizes = sizes[sizes > 0]

         return output
      return None

   def GetParameterSetInfo( self , name ):
      """ dd.shotfile.GetParameterSetInfo( ParameterSetName )\n\nReturns information about the specified parameter set."""
      output = dd_info
      if self.__status:
         info  = self.GetObjectValue( name , 'items' )
         
         error    = ct.c_int32(0)
         _error   = ct.byref( error )
         _diaref  = ct.byref( self.__diaref )
         par_name = ct.c_char_p(name)
         nrec     = ct.c_int32(info)
         _nrec    = ct.byref( nrec )
         rname    = ct.c_char_p(' '*8*info)
         items    = (ct.c_uint32*info)()
         _items   = ct.byref(items)
         format   = (ct.c_uint32*info)()
         _format  = ct.byref( format )
         devsig   = (ct.c_int32*info)()
         _devsig  = ct.byref(devsig)
         lname    = ct.c_ulonglong( len(name) )
         lrname   = ct.c_ulonglong( 8*info )

         result = libddww.ddprinfo_(_error, _diaref, par_name, _nrec, rname, \
                                    _items, _format, _devsig, lname, lrname)
         output.error = error.value
         if error.value != 0:
            self.GetError( error )
         output.N_items = nrec.value
         output.names  = []
         output.items  = []
         output.format = []
         output.devsig = []
         for i in range(info):
            tmp = rname.value[8*i:8*(i+1)]
            if tmp.strip() != '':
               output.names.append(tmp)
         n_pars = len(output.names)
         for j in range(n_pars):
            output.items.append(items[j])
            output.format.append(format[j])
            output.devsig.append(devsig[j])
         
      return output

   def GetParameterInfo(self, set_name, par_name):
      """ Returns information about the parameter 'par_name' of the parameter set 'set_name'."""
      if self.__status:
         output = dd_info
         error   = ct.c_int32(0)
         _error  = ct.byref(error)
         _diaref = ct.byref(self.__diaref)
         pset    = ct.c_char_p(set_name)
         par     = ct.c_char_p(par_name)
         item    = ct.c_uint32(0)
         _item   = ct.byref(item)
         format  = ct.c_uint16(0)
         _format = ct.byref(format)
         lpar = ct.c_ulonglong(len(par_name))
         lset = ct.c_ulonglong(len(set_name))

         result = libddww.dd_prinfo_(_error, _diaref, pset, par, _item, _format, lset, lpar)

         self.GetError(error)
         output.error = error.value
         if error.value == 0:
            output.item = item.value
            output.format = format.value
         return output
      return None

   def GetParameter(self, set_name, par_name):
      """ Returns the value of the parameter 'par_name' of the parameter set 'set_name'. """
      if self.__status:
         info = self.GetParameterInfo(set_name, par_name)
         if info.error == 0:
            error = ct.c_int32(0)
            _error = ct.byref(error)
            _diaref = ct.byref(self.__diaref)
            setn = ct.c_char_p(set_name)
            lset = ct.c_ulonglong(len(set_name))
            par = ct.c_char_p(par_name)
            lpar = ct.c_ulonglong(len(par_name))
            physunit = ct.c_int32(0)
            _physunit = ct.byref(physunit)
# Characters
            if info.format == 2:
               ndim = 1
            if info.format == 1794:
               ndim = 8
            elif info.format == 3842:
               ndim = 16
            elif info.format == 7938:
               ndim = 32
            elif info.format == 16130:
               ndim = 64
            if info.format in (2, 1794, 3842, 7938, 16130):
               nlen = ndim*info.item
               typin = ct.c_int32(6)
               lbuf = ct.c_uint32(nlen)
               buffer = ct.c_char_p('d'*nlen)
               _typin = ct.byref(typin)
               _lbuf = ct.byref(lbuf)
               result = libddww.ddparm_(_error,_diaref,setn,par,_typin,_lbuf,buffer,_physunit,lset,lpar,ct.c_int32(ndim))
               self.GetError(error)
               a=[]
               for j in range(info.item):
                  a.append(buffer.value[j*ndim:(j+1)*ndim])
               return np.array(a)
            else:
               typin = ct.c_int32(fmt2type[info.format])
               lbuf = ct.c_uint32(info.item)
               buffer = (fmt2ct[info.format]*info.item)()
               _typin = ct.byref(typin)
               _lbuf = ct.byref(lbuf)
               _buffer = ct.byref(buffer)
               result = libddww.ddparm_(_error, _diaref, setn, par, _typin, \
                                        _lbuf, _buffer, _physunit, lset, lpar)
               return np.frombuffer(buffer, dtype=np.dtype(buffer))[0]
      return None

   def GetSignal(self, signame, cal=False,nbeg=0,nend=np.infty):
      """Returns the specified signal or sigmal group
      
      :str signame:  the name of the loaded signal
      :bool cal:  True - return calibrated signal (if is avalible, else uncalibrated)
      :int nbeg: return data from the specific index (in case of multidimensional arrays result dependes on the 
                 order of the dimensions
      :int nend: return data to the specific index 

      """
      
      if self.__status:
         info = self.GetInfo(signame)
         self.units = None
         if info.error == 0:
             
            k1 = ct.c_uint32( max(1,nbeg) )
            _k1 = ct.byref( k1 )
            k2 = ct.c_uint32( min(info.ind[0],nend)  )
            _k2 = ct.byref( k2 )
            info.ind[0] = np.subtract(k2,k1)+1
             
            if cal:
               fmt = 2
            else:
               fmt = fmt2type[info.fmt]
            if fmt == 6:
               leng = fmt2len[info.fmt]
            else:
               leng = info.ind[0]

            _diaref  = ct.byref(self.__diaref)
            error    = ct.c_int32(0)
            _error   = ct.byref(error)
            length   = ct.c_uint32(0)
            _length  = ct.byref(length)

            
            
            buffer = self.GetArray(fmt, info.ind)
            _buffer  = ct.byref(buffer)
            cfmt     = ct.c_uint32(fmt)
            _type    = ct.byref(cfmt)
            lbuf = ct.c_uint32(leng)
            _lbuf    = ct.byref(lbuf)
            signam   = ct.c_char_p(signame)
            physdim  = 8*'p'
            _physdim = ct.c_char_p(physdim)
            ncal     = ct.c_int32(0)
            _ncal    = ct.byref(ncal)
            lsig     = ct.c_ulonglong(len(signame))
            lphysdim = ct.c_ulonglong(len(physdim))

            if info.objtyp == 7:
# Signal
               if cal:
# Calibrated Signal
                  result = libddww.ddccsgnl_(_error, _diaref, signam, _k1, _k2, \
                                             _type, _lbuf, _buffer, _length, _ncal, \
                                             _physdim, lsig, lphysdim)
                  self.GetError(error)
                  if error.value != 0:
                     if error.value == 555352323:
                        print 'No Calibration Data, returning uncalibrated data.'
                        return self.GetSignal(signame, cal=False)
                  else:
                     self.units = _physdim.value
               else:
                  result = libddww.ddsignal_(_error, _diaref, signam, _k1, _k2, \
                                             _type, _lbuf, _buffer, _length, lsig)
               return np.frombuffer(buffer, dtype=np.dtype(buffer))[0]

            elif info.objtyp == 6:
# SignalGroup
               if cal:
# Calibrated SignalGroup
                  result = libddww.ddccsgrp_(_error, _diaref, signam, _k1, _k2, \
                                             _type, _lbuf, _buffer, _length, _ncal, \
                                             _physdim, lsig, lphysdim)
                  self.GetError(error)
                  if error.value != 0:
                     if error.value == 556204039:
                        print 'No Calibration Data, returning uncalibrated data.'
                        return self.GetSignal(signame, cal=False)
                  else:
                     self.units = _physdim.value
               else:
                  result = libddww.ddsgroup_(_error, _diaref, signam, _k1, _k2, \
                                             _type, _lbuf, _buffer, _length, lsig)

                  self.GetError(error)
               ind2 = []
               for el in info.ind:
                  if el != 1:
                     ind2.append(el)
               result =  np.reshape(np.frombuffer(buffer, dtype=np.dtype(buffer))[0], newshape=ind2, order='F')
               return result
            else:
               print('The object %s is neither a SIGNAL nor a SIGNAL_GROUP' %signame)
         else:
            print('Error getting SignalGroup %s' %signame)
      return None

   def GetAreabase(self, signame,nbeg=1,nend=np.infty):
      """ 
      Returns the areabase to the specified signal, time is always the first independent variable.
           
      :str signame:  the name of the loaded signal
      :int nbeg: return data from the specific index (in case of multidimensional arrays result dependes on the 
                 order of the dimensions
      :int nend: return data to the specific index 
      
      """
      if self.__status:
         info = self.GetInfo(signame)
         if info.objtyp == 0:
            return None
         if info.error == 0:
            
            k1 = ct.c_uint32( max(1,nbeg) )
            _k1 = ct.byref( k1 )
            k2 = ct.c_uint32( min(info.atlen,nend)  )
            _k2 = ct.byref( k2 )
            
            info.atlen = np.subtract(k2,k1)+1
            #info.sizes[0] = np.subtract(k2,k1)+1
            
            if len(info.sizes)  > 0:
              ind = np.append(info.sizes, info.atlen)
              lbu = info.sizes[0]
            else:
              ind = [info.atlen]
              lbu = 1
                

            

            
            fmt = fmt2type[info.afmt]
            buffer  = self.GetArray(fmt, ind)
            _buffer = ct.byref(buffer)
            error   = ct.c_int32(0)
            _error  = ct.byref(error)
            _diaref = ct.byref(self.__diaref)
            name    = ct.c_char_p(signame)

            lbuf    = ct.c_uint32(lbu)
            _lbuf   = ct.byref(lbuf)
            cfmt    = ct.c_int32(fmt)
            _type   = ct.byref(cfmt)
            leng    = ct.c_int32(0)
            _leng   = ct.byref(leng)
            lname   = ct.c_ulonglong(len(signame))

            result = libddww.ddagroup_(_error, _diaref, name, _k1, _k2, _type, \
                                       _lbuf, _buffer, _leng, lname)

            self.GetError(error)
            ind2 = []
            for el in ind:
               if el != 1 and el != 0:
                  ind2.append(el)
            if len(ind2) > 1:
               abase = np.reshape(np.frombuffer(buffer, dtype = np.dtype(buffer))[0], \
                               newshape=ind2, order='F')
               return abase.T
            elif len(ind2) == 1:
               return np.frombuffer(buffer, dtype = np.dtype(buffer))[0]

      return None
  
  
  
   def GetTimebase(self, signame ,nbeg=1,nend=np.infty):
      """
      Returns the timebase corresponding to the specified signal.
      
      :str signame:  the name of the loaded signal
      :int nbeg: return data from the specific index 
      :int nend: return data to the specific index 

      """
      if self.__status:
         info = self.GetInfo(signame)
         if info.error == 0:
            if info.tlen is None:
               print('No TimeBase found for signal %s' %signame)
               return None
            fmt = fmt2type[info.tfmt]
            
            k1 = ct.c_uint32( max(1,nbeg) )
            _k1 = ct.byref( k1 )
            k2 = ct.c_uint32( min(info.tlen,nend)  )
            _k2 = ct.byref( k2 )
            info.tlen = np.subtract(k2,k1)+1
            buffer = self.GetArray(fmt, [info.tlen])
            _buffer = ct.byref(buffer)
            error = ct.c_int32(0)
            _error = ct.byref(error) 
            _diaref = ct.byref(self.__diaref)
            signam = ct.c_char_p(signame)

            cfmt = ct.c_uint32(fmt)
            _type = ct.byref(cfmt)
            lbuf = ct.c_uint32(info.tlen)
            _lbuf = ct.byref(lbuf)
            length = ct.c_uint32(0)
            _length = ct.byref(length)
            lsig = ct.c_ulonglong(len(signame))

            result = libddww.ddtbase_(_error, _diaref, signam, _k1, _k2, _type, \
                                      _lbuf, _buffer, _length, lsig)
            self.GetError(error)
            if info.tlen > 1:
               return np.frombuffer(buffer, dtype = np.dtype(buffer))[0]
            elif info.tlen == 1:
               return np.frombuffer(buffer, dtype = np.dtype(buffer))
            else:
               print('Array size <= 0')
               return None
      return None

   def GetArray(self, fmt, ind):
      """dd.shotfile.GetArray(fmt, ind, typin=None)\n\nCreates a ct array with the type specified in typ and the format specified in ind."""
      if ind[0] > 0:
         length = 1
         for i in ind:
            if i != 0: length *= i
         if fmt in type2ct.iterkeys():
            return (type2ct[fmt]*int(length))()
         else:
            print('Type %d not supported' %fmt)
      return None

   def GetObjectValue(self, name, field):
      """dd.shotfile.GetObjectValue(Name, Field)\n\nReturns the value specified in Field of the object Name."""
      if self.__status:
         ovalue = ct.c_int32(0)
         if field in ('relations', 'special'):
            ovalue = (ct.c_int32*8)()
         if field in ('format', 'indices'):
            ovalue = (ct.c_int32*3)()
         if field == 'dataformat':
            ovalue = ct.c_uint16(0)
         if field == 'text':
            ovalue = (ct.c_char*64)()
         _value  = ct.byref(ovalue)
         error   = ct.c_int32(0)
         _error  = ct.byref(error)
         _diaref = ct.byref(self.__diaref)
         _name   = ct.c_char_p(name)
         _field  = ct.c_char_p(field)
         lname   = ct.c_ulonglong(len(name))
         lfield   = ct.c_ulonglong(len(field))

         result = libddww.ddobjval_(_error,_diaref,_name,_field,_value,lname,lfield)

         self.GetError(error)
         if np.size(ovalue) == 1:
            return ovalue.value
         else:
            return np.frombuffer(ovalue, dtype = np.dtype(ovalue))[0]
      return None

   def GetObjectData(self, obj_name):
      """dd.shotfile.GetObjectData(ObjectName, fmt=ct.c_float)\n\nReturns the data part of the specified object and interpretes it as an array of the datatype specified in fmt."""
      if self.__status:
         length = self.GetObjectValue(obj_name, 'length') 
         nsize  = self.GetObjectValue(obj_name, 'size')
         fmt    = self.GetObjectValue(obj_name, 'dataformat')
         cfmt = fmt2ct[fmt]

         error   = ct.c_int32(0)
         _error  = ct.byref(error)
         _diaref = ct.byref(self.__diaref)
         name    = ct.c_char_p(obj_name)
         clength = ct.c_int32(length)
         _length = ct.byref(clength)
         buffer  = (ct.c_byte*length)()
         _buffer = ct.byref(buffer)
         length  = ct.c_int32(0)
         _length = ct.byref(leng)
         lname   = ct.c_ulonglong(len(obj_name))

         result = libddww.ddobjdata_(_error, _diaref, name, _length, _buffer, \
                                     _length, lname)

         self.GetError(error)
         darr = np.frombuffer((cfmt*(length/ct.sizeof(cfmt))).from_buffer(buffer), dtype = np.dtype(cfmt))
         ny = len(darr)/nsize
         if ny == 1:
            return darr
         else:
            return darr.reshape(nsize, ny, order='F')
      return None

   def GetListInfo(self):
      if self.__status:
         error = ct.c_int32(0)
         _error = ct.byref(error)
         _diaref = ct.byref(self.__diaref)
         lbuf = ct.c_uint32(255)
         _lbuf = ct.byref(lbuf)
         buf = (ct.c_int32*255)()
         _buf = ct.byref(buf)
         length = ct.c_int32(0)
         _length = ct.byref(length) 
         result = libddww.ddflist_(_error, _diaref, _lbuf, _buf, _length)
         self.GetError(error)
         if length.value != 0:
            return np.int(buf[0:np.int(length)])
      return None

   def GetObjectHeader(self, name):
      if self.__status:
         output = dd_info
         error = ct.c_int32(0)
         _error = ct.byref(error)
         _diaref = ct.byref(self.__diaref)
         _name = ct.c_char_p(name)
         lname = ct.c_ulonglong(len(name))
         buffer = (ct.c_int32*26)()
         _buffer = ct.byref(buffer)
         text = 64*'t'
         _text = ct.c_char_p(text)
         ltext = ct.c_ulonglong(len(text))
         result = libddww.ddobjhdr_(_error, _diaref, _name, _buffer, _text, lname, ltext)
         self.GetError(error)
         output.error = error.value
         if error.value == 0:
            output.buffer = buffer[0:26]
            output.text   = _text.value
         return output
      return None

   def GetObjectName(self, obj):
      if self.__status:
         name = 8*'1'
         error   = ct.c_int32(0)
         _error  = ct.byref(error)
         _diaref = ct.byref(self.__diaref)
         _name   = ct.c_char_p(name)
         obje    = ct.c_int32(obj)
         _obje   = ct.byref(obje)
         lname   = ct.c_ulonglong(len(name))

         result = libddww.ddobjname_(_error, _diaref, _obje, _name, lname)

         if error.value == 0:
            return _name.value.replace(' ','')
      return -1

   def GetNames(self):
      if self.__status:
         i = 1
         ok = True
         result = []
         while self.GetObjectName(i) != -1:
            result.append(self.GetObjectName(i))
            i += 1
         if i != 0:
            return result
      return None

   def GetDescr(self,name):
      descr = ''
      for str in self.GetObjectValue(name,'text'):
         descr += str
      return descr

   def GetRelations(self, name):

      if self.__status:
         rel_out = dd_info
         rel_out.error = 0
         ids = self.GetObjectValue(name, 'relations')
         if np.max(np.abs(ids)) == 0:
            rel_out.error = 1
         else:
            rel_out.id  = []
            rel_out.typ = []
            rel_out.txt = []
            for objid in ids:
               if objid != 65535:
                  rel_out.id.append(objid)
                  tname = self.GetObjectName(objid)
                  rel_out.typ.append(self.GetObjectValue(tname, 'objtype'))
                  rel_out.txt.append(tname)
         return rel_out
      return None

   def cShotNr(self, experiment, diagname, shotnumber):
      """ dd.shotfile.cShotNr(experiment, diagname, shotnumber) """
      exp = ct.c_char_p(experiment)
      diag = ct.c_char_p(diagname)
      cshot = ct.c_uint32(0)
      _cshot = ct.byref(cshot)
      shot = ct.c_uint32(shotnumber)
      _shot = ct.byref(shot)
      lexp = ct.c_ulonglong(len(experiment))
      ldiag = ct.c_ulonglong(len(diagname))
      result = libddww.ddcshotnr_(exp,diag, _shot, _cshot, lexp, ldiag)
      self.GetError(ct.c_int32(result))
      return cshot.value

   def LastShotNr(self):
      """ nshot = dd.shotfile.LastShotNr() """
      error  = ct.c_int(0)
      _error = ct.byref(error)
      cshot  = ct.c_uint32(0)
      _cshot = ct.byref(cshot)
      result = libddww.ddlastshotnr_(_error, _cshot)
      return cshot.value

   def GetObject(self, name, cal=False,nbeg=0,nend=np.infty):

      info = self.GetInfo(name)

      output = dd_info
      output.data  = None
      output.area  = None
      output.time  = None
      output.units = self.units
      output.tname = info.tname
      output.aname = info.aname
      
      if output.objtyp not in (6, 7, 8, 13):
         print('Only TB, AB, SIG or SGR supported')
         return None

      if output.tname != None:
         output.time  = self.GetTimebase(name,nbeg=nbeg,nend=nend)

      if output.objtyp == 8:
         output.data  = self.GetTimebase(name,nbeg=nbeg,nend=nend)

      if output.objtyp == 13:
         output.data  = self.GetAreabase(name,nbeg=nbeg,nend=nend)

      if output.objtyp == 7:
         output.data  = self.GetSignal(name, cal=cal,nbeg=nbeg,nend=nend)

      if output.objtyp == 6:
         output.data  = self.GetSignal(name, cal=cal,nbeg=nbeg,nend=nend)
         if output.aname != None:
            output.area  = self.GetAreabase(name)
            if output.index != None:
               if output.index in (2, 3):
                  output.data = output.data.T
         if output.area != None:
            arel = self.GetRelations(info.aname)
            if len(output.area.shape) == 1:
               if 8 not in arel.typ:
                  output.area = np.tile(output.area, (len(output.time),1))
               else:
                  output.area = np.reshape(output.area, newshape=(len(output.area),1))

      return output

# Backward compatibility

   def GetSignalCalibrated(self, name,nbeg=1,nend=np.infty):
      """ 
      Return  calibrated 1D signal
           
      :str signame:  the name of the loaded signal
      :int nbeg: return data from the specific index
      :int nend: return data to the specific index 
      
      """
      return self.GetSignal(name, cal=True,nbeg=nbeg,nend=nend)

   def GetSignalGroup(self, name,nbeg=1,nend=np.infty):
      """ 
       Return  uncalibrated nD signal
           
      :str signame:  the name of the loaded signal
      :int nbeg: return data from the specific index (in case of multidimensional arrays result dependes on the 
                 order of the dimensions
      :int nend: return data to the specific index 
      
      """
      return self.GetSignal(name,nbeg=nbeg,nend=nend)

   def GetSignalGroupCalibrated(self, name,nbeg=1,nend=np.infty):
      """ 
       Return  calibrated nD signal

           
      :str signame:  the name of the loaded signal
      :int nbeg: return data from the specific index (in case of multidimensional arrays result dependes on the 
                 order of the dimensions
      :int nend: return data to the specific index 
      
      """
      return self.GetSignal(name, cal=True,nbeg=nbeg,nend=nend)

   def GetSignalInfo(self,name):
      return self.GetInfo(name)
   def GetTimebaseInfo(self,name):
      return self.GetInfo(name)
  
  
class Ds_help:
   status=False

class Dataset:

   def __init__(self, diag, nshot, exp='AUGD', ed=0, cal=False):

      self.variables = {}
      sf = shotfile()

      if sf.Open(diag, nshot, experiment=exp, edition=ed):

         lis = sf.GetNames()

         for obj in lis:
            self.variables[obj] = Ds_help()
            info = sf.GetInfo(obj)
            self.variables[obj].devtyp = obj_d[info.objtyp]
            if info.objtyp in (6, 7, 8, 13):
               tmp = sf.GetObject(obj, cal=cal)
               self.variables[obj].data = tmp.data
               self.variables[obj].dimensions = []
               if tmp.tname != None:
                   self.variables[obj].dimensions.append(tmp.tname)
               if tmp.aname != None:
                   self.variables[obj].dimensions.append(tmp.aname)
               self.variables[obj].units = tmp.units
               self.variables[obj].long_name = sf.GetDescr(obj)
               rels = sf.GetRelations(obj).txt
               self.variables[obj].relations = rels 
            if info.objtyp == 4:
               ps_info = sf.GetParameterSetInfo(obj)
               self.variables[obj].data = {}
               for pn in ps_info.names: 
                  self.variables[obj].data[pn] = sf.GetParameter(obj,pn)
         sf.Close()
