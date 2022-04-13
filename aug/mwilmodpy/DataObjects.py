# MDSplus based data objects

import pmds as mds
import numpy as np


# MDSplus data source:
#   server and tree (common for many signals)

class DataSource:

  def __init__ (self, server, treename):
    global lastserver, lasttree, lastshot
    self.server = server
    self.treename = treename

# (keep) open connection/tree for shot and return True if successful
  def open(self, shot):
    global lastserver, lasttree, lastshot
    if self.server<>lastserver:
      print 'Connecting to MDSplus server '+self.server
      try:
        mds.mdsconnect(self.server)
        lastserver = self.server
      except RuntimeError:
        print "ERROR when trying to connect to: ", server
        lastserver = ''
        return False
    if self.treename<>'':
      if (self.treename<>lasttree) | (shot<>lastshot):
        print 'Open tree '+self.treename+' for shot ', shot
        try:
          mds.mdsvalue('treeopen("'+self.treename+'",'+str(shot)+')')
          lasttree = self.treename
          lastshot = shot
        except RuntimeError:
          print "ERROR opening tree ", self.treename, " for shot ", shot
          lasttree = ''
          lastshot = -1
          return False
    return True

# ------------------------------------------------------------


# anything that can be read with a single TDI expression

class Expression:

  def __init__ (self, source, name, expression):
    self.source = source
    self.lastshot = -1
    self.havedata = False
    self.name = name
    self.expression = expression
    
  def read(self, shot):
    if self.source.open(shot):
      try:
        self.data = mds.mdsvalue(self.expression)
      except RuntimeError:
        print "ERROR reading data of ", self.expression
        self.havedata = False
        return False
      self.lastshot = shot
      self.havedata = True
      return True


# generic Time trace

class TimeTrace:

  def __init__ (self, source, name, expression, fallback1=None, fallback2=None, channel=None):
    self.source = source
    self.lastshot = -1
    self.havedata = False
    self.name = name
    self.expression = expression
    self.fallback1 = fallback1
    self.fallback2 = fallback2
    self.channel = channel

  def read(self, shot):
    if self.source.open(shot):
      try:
        data = mds.mdsvalue('_s=('+self.expression+')')
      except RuntimeError:
        print "ERROR reading data of ", self.expression

        if self.fallback1 is None:
          self.havedata = False
          return False
        print "Falling back to ", self.fallback1
        try:
          data = mds.mdsvalue('_s=('+self.fallback1+')')
        except RuntimeError:
          print "ERROR reading data of ", self.fallback1

          if self.fallback2 is None:
            self.havedata = False
            return False
          print "Falling back to ", self.fallback2
          try:
            data = mds.mdsvalue('_s=('+self.fallback2+')')
          except RuntimeError:
            print "ERROR reading data of ", self.fallback2
            self.havedata = False
            return False

      if self.channel is None:
        self.data = data
      else:
        self.data = data[self.channel-1,:]

      try:
        self.timebase = mds.mdsvalue('dim_of(_s)')
      except RuntimeError:
        print "ERROR reading dimension of ", self.expression
        self.havedata = False
        return False
      self.lastshot = shot
      self.havedata = True
      return True



# profile vs. R,z,t
# values are from expression 'expr' which is always evaluated first
# R, z, and time are from 'Rexpr', 'Zexpr', 'Texpr' (in this order)

class ProfileRzt:

  def __init__ (self, source, name, expr, Rexpr, Zexpr, Texpr):
    self.source = source
    self.lastshot = -1
    self.havedata = False
    self.name = name
    self.expression = expr
    self.Rexpr = Rexpr
    self.Zexpr = Zexpr
    self.Texpr = Texpr

  def read(self, shot):
    if self.source.open(shot):
      try:
        self.data = mds.mdsvalue(self.expression)
      except RuntimeError:
        print "ERROR reading data of ", self.expression
        self.havedata = False
        return False
      try:
        self.R = mds.mdsvalue(self.Rexpr)
      except RuntimeError:
        print "ERROR reading R coordinate through ", self.Rexpr
        self.havedata = False
        return False
      try:
        self.Z = mds.mdsvalue(self.Zexpr)
      except RuntimeError:
        print "ERROR reading Z coordinate through ", self.Zexpr
        self.havedata = False
        return False
      try:
        self.timebase = mds.mdsvalue(self.Texpr)
      except RuntimeError:
        print "ERROR reading time base through ", self.Texpr
        self.havedata = False
        return False
      self.lastshot = shot
      self.havedata = True
      return True


# profile vs. Psi
# values are from expression 'expr' which is always evaluated first
# Psi and time are from 'Psiexpr', 'Texpr' (in this order)

class ProfilePsi:

  def __init__ (self, source, name, Dataexpr, Errexpr, Psiexpr, Texpr):
    self.source = source
    self.lastshot = -1
    self.havedata = False
    self.name = name
    self.expression = Dataexpr
    self.Errexpr = Errexpr
    self.Psiexpr = Psiexpr
    self.Texpr = Texpr

  def read(self, shot):
    if self.source.open(shot):

      if self.expression<>"":
        try:
          self.data = mds.mdsvalue(self.expression)
        except RuntimeError:
          print "ERROR reading data of ", self.expression
          self.havedata = False
          return False

      if self.Errexpr<>"":
        try:
          self.errorbar = mds.mdsvalue(self.Errexpr)
        except RuntimeError:
          print "ERROR reading data of ", self.Errexpr
          self.havedata = False
          return False

      if self.Psiexpr<>"":
        try:
          self.Psi = mds.mdsvalue(self.Psiexpr)
        except RuntimeError:
          print "ERROR reading Psi coordinate through ", self.Psiexpr
          self.havedata = False
          return False

      if self.Texpr<>"":
        try:
          self.timebase = mds.mdsvalue(self.Texpr)
        except RuntimeError:
          print "ERROR reading time base through ", self.Texpr
          self.havedata = False
          return False

      self.lastshot = shot
      self.havedata = True
      return True


# initialise this module

lastserver = ''
lasttree = ''
lastshot = -1
