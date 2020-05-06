import ctypes as ct
import numpy as np

unit_d = { \
  0:None, 1:'kg', 2:'m', 3:'V', 4:'A', 5:'mV', 6:'eV', 7:'J', 8:'s', 9:'min', 10:'h', \
  11:'Celsius', 12:'pm', 13:'msec', 14:'1/V', 15:'K', 16:'degree', 17:'keV', 18:'cm', \
  19:'mm', 20:'micron', 21:'+-5V/12b', 22:'+-10V/12', 23:'counts', 24:'10e14/cc', \
  25:'Vs', 26:'A/(m*m)', 27:'T', 28:'W', 29:'C', 30:'m^2', 31:'m^3', 32:'kA', \
  33:'W/m^2', 34:'W/m^2/nm', 35:'1/m', 36:'1/m^2', 37:'1/m^3', 38:'10e19/m^', \
  39:'mbar', 40:'Pa', 41:'bar', 42:'kV', 43:'mA', 44:'+-5V/409', 45:'+-10V/40', \
  46:'Hz', 47:'+5V/4095', 48:'+10V/409', 49:'l/min', 50:'1/s', 51:'MN/m', 52:'MJ', \
  53:'ASCII', 54:'V/A', 55:'m^3/h', 56:'MW', 57:'mm^2/s', 58:'m^2/s', 59:'W/(mm*K)', \
  60:'1/mm', 61:'dB', 62:'1/J', 63:'MW/m^2', 64:'kW/m^2', 65:'kA/s', 66:'T/s', \
  67:'W/(m^2*s', 68:'W/m^3', 69:'cnts/s', 70:'m/s', 71:'rad/s', 72:'GHz', 73:'N/A', \
  74:'nm', 75:'+-5V/16b', 76:'+-10V/16', 77:'AU', 78:'kW', 79:'J/m^2', 80:'V/m', \
  81:'Ph/(qm*s', 82:'1/(m^2*s', 83:'kA^2*s', 84:'Nm', 85:'+5V/12bi', 86:'+10V/12b', \
  87:'+-5V/13b', 88:'+-10V/13', 89:'+5V/13bi', 90:'+10V/13b', 91:'+-5V/819', \
  92:'+-10V/81', 93:'+5V/8191', 94:'+10V/819', 95:'+-5V/14b', 96:'+-10V/14', \
  97:'+5V/14bi', 98:'+10V/14b', 99:'+-5V/163', 100:'+-10V/16', 101:'+5V/1638', \
  102:'+10V/163', 103:'+-5V/15b', 104:'+-10V/15', 105:'+5V/15bi', 106:'+10V/15b', \
  107:'+-5V/327', 108:'+-10V/32', 109:'+5V/3276', 110:'+10V/327', 111:'+5V/16bi', \
  112:'+10V/16b', 113:'+-5V/655', 114:'+-10V/65', 115:'+5V/6553', 116:'+10V/655', \
  117:'nanosec', 118:'amu', 119:'pct', 120:'MHz' \
}

obj_d = { 1: 'Diagnostic', 2: 'List',       3: 'Device',      4: 'Param_Set', \
          5: 'Map_Func',   6: 'Sig_Group',  7: 'Signal',      8: 'Time_Base', \
          9: 'SF_List',    10: 'Algorithm', 11: 'Update_Set', 12: 'Loc_Timer', \
          13: 'Area_Base', 14: 'Qualifier', 15: 'ModObj',     16: 'Map_Extd',  \
          17: 'Resource'}

# Source: /afs/ipp/aug/ads/common/codes/formats
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

fmt2type = { 2:6, 3:11, 4:1, 5:2, 6:3, 7:5, 9:12, 13:10, 14:13, 15:14, \
             1794:6, 3842:6, 7938:6, 12034:6, 16130:6, 18178:6}
fmt2len = {2:1, 1794:8, 3842:16, 7938:32, 12034:48, 16130:64, 18178:72}

# Corresponding ctypes datatypes
# FIXME: Raw (0), Complex (4)

type2ct = { 1:ct.c_int32, 2:ct.c_float, 3:ct.c_double, 5: ct.c_bool, 6:ct.c_char, \
  10: ct.c_int64, 11:ct.c_int16, 12:ct.c_uint16, 13:ct.c_uint32, 14:ct.c_uint64 }

# FIXME: Raw (0), Complex (4)

type2np = { 1:np.int32, 2:np.float, 3:np.double, 6:np.char, \
  10: np.int64, 11:np.int16, 12:np.uint16, 13:np.uint32, 14:np.uint64 }

fmt2ct = { 2:ct.c_char, 3:ct.c_int16, 4:ct.c_int32, 5:ct.c_float, \
   6:ct.c_double, 9:ct.c_uint16, 13:ct.c_int64, 14:ct.c_uint32, 15:ct.c_uint64 }

fmt2np = { 2:np.char, 3:np.int16, 4:np.int32, 5:np.float32, \
   6:np.float64, 9:np.uint16, 13:np.int64, 14:np.uint32, 15:np.uint64 }
