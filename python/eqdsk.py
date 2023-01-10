#!/usr/bin/env python3
# -*- coding: utf-8 -*-

class eqdsk_file:

  def __init__ (self, filename:str):
    self.generate_from_file(filename)

    self.verbose = False


  def generate_from_file(self, filename:str):
    """
    \note Might not be suitable for general format of eqdsk file.
    """
    from numpy import array

    with open(filename) as f:
      l = f.readlines()

    l = [ls.strip() for ls in l]

    self.filenamedateEQDSK = l[0][0:48]
    # l[0].split()[1] is dummy, not used
    # number r coordinates for the grid
    self.nrgr = int(l[0].split()[2]);
    # number z coordinates for the grid
    self.nzgr = int(l[0].split()[3]);
    # Read boxlength
    self.rboxlength = float(l[1].split()[0])
    self.zboxlength = float(l[1].split()[1])
    self.R0 = float(l[1].split()[2])
    self.rboxleft = float(l[1].split()[3])
    self.zbox_mid = float(l[1].split()[4])
    self.R = [self.rboxleft + self.rboxlength/(self.nrgr-1)*k for k in range(self.nrgr)]
    self.Z = [self.zbox_mid - self.zboxlength/2 + self.zboxlength/(self.nzgr-1)*k for k in range(self.nzgr)] # Needs to be transposed.
    # Read magnetic axis
    self.Rpsi0 = float(l[2].split()[0])
    self.Zpsi0 = float(l[2].split()[1])
    # Read psiaxis and psiedge
    self.PsiaxisVs = float(l[2].split()[2])
    self.PsiedgeVs = float(l[2].split()[3])
    # Read B0 and current
    self.Btor_at_R0=float(l[2].split()[4])
    self.Ip=float(l[3].split()[0])
    # Rest of line 3 (4 one-based) and the next line is ignored.

    startline = 5 # 0-based index
    # Read T(psi) - f (size nrgr)
    self.psieqprof = [1.0/(self.nrgr-1)*k for k in range(self.nrgr)]
    self.fprof = []
    line_number = startline
    while (len(self.fprof) < self.nrgr):
      self.fprof += [float(n) for n in l[line_number].split()]
      line_number += 1

    # Read p(psi) (size nrgr)
    self.ptotprof = []
    if self.nrgr % 5 == 0:
      line_number += 1
    while (len(self.ptotprof) < self.nrgr):
      self.ptotprof += [float(n) for n in l[line_number].split()]
      line_number += 1

    # Read TT'(psi) - ff' (size nrgr)
    self.fdfdpsiprof = []
    if self.nrgr % 5 == 0:
      line_number += 1
    while (len(self.fdfdpsiprof) < self.nrgr):
      self.fdfdpsiprof += [float(n) for n in l[line_number].split()]
      line_number += 1

    # Read p'  (size nrgr)
    self.dpressdpsiprof = []
    if self.nrgr % 5 == 0:
      line_number += 1
    while (len(self.dpressdpsiprof) < self.nrgr):
      self.dpressdpsiprof += [float(n) for n in l[line_number].split()]
      line_number += 1

    # Read psi (size nrgr x nzgr in file?, nzgr x nrgr in 'output')
    self.PsiVs = []
    if self.nrgr % 5 == 0:
      line_number += 1
    while (len(self.PsiVs) < self.nrgr*self.nzgr):
      self.PsiVs += [float(n) for n in l[line_number].split()]
      line_number += 1
    self.psinorm = [(PsiV-self.PsiaxisVs+self.PsiedgeVs)/(-self.PsiaxisVs+self.PsiedgeVs) for PsiV in self.PsiVs]
    self.PsiVs = array(self.PsiVs)
    self.PsiVs = self.PsiVs.reshape(self.nzgr,self.nrgr)

    # Read Q profile
    self.qprof = []
    if (self.nrgr*self.nzgr) % 5 == 0:
      line_number += 1
    while (len(self.qprof) < self.nrgr):
      self.qprof += [float(n) for n in l[line_number].split()]
      line_number += 1

    # Read npbound and nblimiter (in one line)
    if self.nrgr % 5 == 0:
      line_number += 1
    self.npbound = int(l[line_number].split()[0])
    self.nplimiter = int(l[line_number].split()[1])
    # Read boundary (LCFS)
    line_number += 1
    self.Lcfs = []
    while (len(self.Lcfs) < 2*self.npbound):
      self.Lcfs += [float(n) for n in l[line_number].split()]
      line_number += 1
    self.Lcfs = array(self.Lcfs)
    self.Lcfs = self.Lcfs.reshape(self.npbound, 2)
    self.Lcfs = self.Lcfs.transpose()

    # Read limiter (LCFS)
    if (2*self.npbound) % 5 == 0:
      line_number += 1
    self.Limiter = []
    while (len(self.Limiter) < 2*self.nplimiter):
      self.Limiter += [float(n) for n in l[line_number].split()]
      line_number += 1
    self.Limiter = array(self.Limiter)
    self.Limiter = self.Limiter.reshape(self.nplimiter, 2)
    self.Limiter = self.Limiter.transpose()

    if (2*self.nplimiter) % 5 == 0:
      line_number += 1
    self.ncoils = int(l[line_number])
    line_number += 1
    # Read PFcoilData
    self.PFcoilData = []
    while (len(self.PFcoilData) < 5*self.ncoils):
      self.PFcoilData += [float(n) for n in l[line_number].split()]
      line_number += 1
    self.PFcoilData = array(self.PFcoilData)
    self.PFcoilData = self.PFcoilData.reshape(self.ncoils, 5)

    self.Rcond = self.PFcoilData[:,0]
    self.Zcond = self.PFcoilData[:,1]
    self.DRcond = self.PFcoilData[:,2]
    self.DZcond = self.PFcoilData[:,3]
    self.Icond = self.PFcoilData[:,4]

    if (5*self.ncoils) % 5 == 0:
      line_number += 1
    # Read Brn
    self.Brn = []
    while (len(self.Brn) < self.nrgr*self.nzgr):
      self.Brn += [float(n) for n in l[line_number].split()]
      line_number += 1
    self.Brn = array(self.Brn)
    self.Brn = self.Brn.reshape(self.nzgr, self.nrgr)

    if (self.nrgr*self.nzgr) % 5 == 0:
      line_number += 1
    # Read Bzn
    self.Bzn = []
    while (len(self.Bzn) < self.nrgr*self.nzgr):
      self.Bzn += [float(n) for n in l[line_number].split()]
      line_number += 1
    self.Bzn = array(self.Bzn)
    self.Bzn = self.Bzn.reshape(self.nzgr, self.nrgr)

    self.rest = l[line_number:]


  def sort_into_coilgroups(self, indices_coilgroupstarts:list, coilgroups:list, coiltags:list):
    """
    This method takes the coils and sorts it into groups according to provided information.

    After reading coil data is just a list of points with currents. This
    method will sort them into coil(-groups).
    As it is not clear if this process can be automated at all, input is
    required to do so.

    input:
    ------
    indices_coilgroupstarts: list of integers, indices of the points
      which represent the start of each coil.
    coilgroups: list of integers, same length as indices_coilgroupstarts.
      An 'ID' for the corresponding coil. Coils with same ID can not be
      changed independent(?).
    coiltags: list of strings, same length as indices_coilgroupstarts.
      Short(?) descriptive name for coilgroup.
    """
    self.coilgroups = []

    for k in range(len(indices_coilgroupstarts)):
      i = indices_coilgroupstarts[k]
      imax = len(self.Rcond)
      if k < len(indices_coilgroupstarts)-1:
        imax = indices_coilgroupstarts[k+1]

      for j in range(i, imax):
        c = Coilgroup(self.Icond[j], coilgroups[k], coiltags[k])

        if c.Tag in {'CS','PF'}:
          p = Point(self.Rcond[j],
                    -self.Rcond[j]/50.0,
                    self.Zcond[j])
          c.coordinates.append(p)

          p = Point(self.Rcond[j],
                    +self.Rcond[j]/50.0,
                    self.Zcond[j])
          c.coordinates.append(p)
        else:
          p = Point(self.Rcond[j] - self.DRcond[j],
                    0.0,
                    self.Zcond[j] - self.DZcond[j])
          c.coordinates.append(p)

          p = Point(self.Rcond[j] - self.DRcond[j],
                    0.0,
                    self.Zcond[j] + self.DZcond[j])
          c.coordinates.append(p)

          p = Point(self.Rcond[j] + self.DRcond[j],
                    0.0,
                    self.Zcond[j] + self.DZcond[j])
          c.coordinates.append(p)

          p = Point(self.Rcond[j] + self.DRcond[j],
                    0.0,
                    self.Zcond[j] - self.DZcond[j])
          c.coordinates.append(p)

        self.coilgroups.append(c)


  def increase_periodicity_in_coilgroups(self, periods:int):
    """
    Increase the number of periods for the coilgroups.

    Increase the periodicity, by increasing the number of elements, i.e.
    the makegrid file will still have 'periods 1'.

    Existing data needs to be rotated. Angle is

      k * 2*pi/periods

    with $k \in [0, periods-1]$ or $k \in [1, periods]$. The case
    $k = 0$ and $k = periods$ corespond to no rotation/one full rotation
    and thus can be ignored.
    The rotation itself is along the z-axis and can thus be described by

      x' = x \cos \alpha  - y \sin \alpha
      y' = x \sin \alpha  + y \cos \alpha
      z' = z

    input:
    ------
    periods: integer, number by which to mulitply number of elements.

    output:
    -------
    None

    sideeffects:
    ------------
    Adds elements to member coilgroups.
    """
    from math import cos, pi, sin

    for c in self.coilgroups.copy():
      coordinates = c.coordinates.copy()
      for k in range(1,periods):
        alpha = float(k) * 2.0*pi/float(periods)

        # Differentiate between elements that have toroidal symmetry and
        # those that have no (continous) toroidal symmetry.
        if c.Tag in {'CS','PF'}:
          for p in coordinates:
            c.coordinates.append(Point(p.x*cos(alpha) - p.y*sin(alpha),
                                       p.x*sin(alpha) + p.y*cos(alpha),
                                       p.z))
        else:
          cn = Coilgroup(c.I, c.coilgroup, c.Tag)
          for p in coordinates:
            cn.coordinates.append(Point(p.x*cos(alpha) - p.y*sin(alpha),
                                        p.x*sin(alpha) + p.y*cos(alpha),
                                        p.z))

          self.coilgroups.append(cn)


  def to_makegrid_input(self, filename):
    """
    Creates a makegrid coil input file from the class.

    This method expects, that the coilgroups are set.

    From makegrid documentation (https://princetonuniversity.github.io/STELLOPT/MAKEGRID):

    The coils file is a standard filament format which has it's origin
    in coil design software. The first line indicates the number of
    field periods for a run. The next line indicates the beginning of
    the filament section. The third line is read by MAKEGRID but
    ignored. Then the definition of a coil begins. The coils is
    specified by connecting points in space along a line defining a
    coil. For each point a current is specified. The last point should
    correspond with the first in a given loop and have zero current. In
    the following example, a square coil is specified.

    periods 3
    begin filament
    mirror NIL
    0.0   0.0   0.0   5.0
    1.0   0.0   0.0   5.0
    1.0   1.0   0.0   5.0
    0.0   1.0   0.0   5.0
    0.0   0.0   0.0   0.0   1   SQ
    end

    The last line also contains a number which associates this coil with
    coilgroup 1 and a name 'SQ.' The fields from all coils of the same
    coilgroup will be summed. In a sense, the user can think of
    coilgroups as indicating common (or coupled) power supplies.

    \todo coordinate system?

    input:
    ------
    filename: string, path+name of the file to which to write the result.

    output:
    -------
    none

    sideeffects:
    ------------
    Creates file on disk and writes results to the file.
    """

    periods = 1

    with open(filename, 'w') as f:
      f.write("periods {}\n".format(periods))
      f.write("begin filament\n")
      f.write("mirror NIL\n")

      for c in self.coilgroups:
        if abs(c.I) < 1.0e-10:
          if self.verbose:
            print("Ignoring coilgroup with current {} and tag {}".format(c.I, c.Tag))
          continue

        for p in c.coordinates:
          f.write("{} {} {} {}\n".format(p.x, p.y, p.z, c.I))

        p = c.coordinates[0]
        f.write("{} {} {} {} {} {}\n".format(p.x, p.y, p.z, 0.0, c.coilgroup, c.Tag))

      f.write("end\n")

class Coilgroup:
  """
  Collect data for a coil(group):
  - current
  - id
  - Tag
  - coordinates
  """
  def __init__(self, I, coilgroup, Tag):
    self.I = I
    self.coilgroup = coilgroup
    self.Tag = Tag
    self.coordinates = []

class Point:
  """
  Collect data for a 3D point, i.e. x, y and z position (not necessarily
  cartesian).
  """
  def __init__(self, x, y, z):
    self.x = x
    self.y = y
    self.z = z
