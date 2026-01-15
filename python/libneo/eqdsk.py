class eqdsk_file:

  def __init__(self, filename: str):
    self.verbose = False
    self._psi_interpolator = None
    self.generate_from_file(filename)

  def _build_psi_interpolator(self):
    from scipy.interpolate import RegularGridInterpolator
    import numpy as np
    self._psi_interpolator = RegularGridInterpolator(
        (self.Z, self.R), self.PsiVs, method='cubic', bounds_error=False
    )
    i_r_axis = np.argmin(np.abs(self.R - self.Rpsi0))
    i_z_axis = np.argmin(np.abs(self.Z - self.Zpsi0))
    self._psi_at_axis = self.PsiVs[i_z_axis, i_r_axis]
    self._psi_at_edge = self.PsiedgeVs
    if np.sign(self._psi_at_axis) != np.sign(self.PsiaxisVs):
      self._psi_at_edge = -self.PsiedgeVs

  def psi_at_rz(self, R, Z, grid=False):
    """
    Interpolate poloidal flux psi at arbitrary (R, Z) coordinates.

    Parameters
    ----------
    R : float or array_like
        Major radius coordinate(s) in meters.
    Z : float or array_like
        Vertical coordinate(s) in meters.
    grid : bool, optional
        If True, evaluate on the grid formed by R and Z arrays.
        If False (default), evaluate at corresponding pairs (R[i], Z[i]).

    Returns
    -------
    psi : float or ndarray
        Poloidal flux value(s) at the given coordinates.
    """
    import numpy as np
    if self._psi_interpolator is None:
      self._build_psi_interpolator()

    R = np.atleast_1d(R)
    Z = np.atleast_1d(Z)

    if grid:
      ZZ, RR = np.meshgrid(Z, R, indexing='ij')
      points = np.column_stack([ZZ.ravel(), RR.ravel()])
      result = self._psi_interpolator(points).reshape(ZZ.shape)
    else:
      points = np.column_stack([Z, R])
      result = self._psi_interpolator(points)
      if result.size == 1:
        result = result.item()

    return result

  def spol_at_rz(self, R, Z):
    """
    Compute normalized poloidal flux coordinate s_pol at (R, Z).

    s_pol = (psi - psi_axis) / (psi_edge - psi_axis)

    Parameters
    ----------
    R : float or array_like
        Major radius coordinate(s) in meters.
    Z : float or array_like
        Vertical coordinate(s) in meters.

    Returns
    -------
    s_pol : float or ndarray
        Normalized poloidal flux, 0 at axis, 1 at separatrix.
    """
    if self._psi_interpolator is None:
      self._build_psi_interpolator()
    psi = self.psi_at_rz(R, Z)
    return (psi - self._psi_at_axis) / (self._psi_at_edge - self._psi_at_axis)

  def theta_geometric_at_rz(self, R, Z):
    """
    Compute geometric poloidal angle at (R, Z) relative to magnetic axis.

    theta = atan2(Z - Z_axis, R - R_axis)

    Parameters
    ----------
    R : float or array_like
        Major radius coordinate(s) in meters.
    Z : float or array_like
        Vertical coordinate(s) in meters.

    Returns
    -------
    theta : float or ndarray
        Geometric poloidal angle in radians, range (-pi, pi].
        theta=0 at outboard midplane, theta=pi/2 at top.
    """
    import numpy as np
    R = np.atleast_1d(R)
    Z = np.atleast_1d(Z)
    theta = np.arctan2(Z - self.Zpsi0, R - self.Rpsi0)
    if theta.size == 1:
      return theta.item()
    return theta

  def rz_to_flux_coords(self, R, Z):
    """
    Convert (R, Z) cylindrical coordinates to flux coordinates (s_pol, theta).

    Parameters
    ----------
    R : float or array_like
        Major radius coordinate(s) in meters.
    Z : float or array_like
        Vertical coordinate(s) in meters.

    Returns
    -------
    s_pol : float or ndarray
        Normalized poloidal flux coordinate.
    theta : float or ndarray
        Geometric poloidal angle in radians.
    """
    s_pol = self.spol_at_rz(R, Z)
    theta = self.theta_geometric_at_rz(R, Z)
    return s_pol, theta


  def generate_from_file(self, filename:str):
    """
    \note Might not be suitable for general format of eqdsk file.
    """
    import numpy as np

    with open(filename, "rb") as f:
      # First line
      self.filenamedateEQDSK = f.read(48).decode("utf-8")
      line = f.readline().split()

      # number of r coordinates for the grid
      self.nrgr = int(line[1])
      # number of z coordinates for the grid
      self.nzgr = int(line[2])

      # Second line
      rewind = 0
      try:
        rewind += 16
        self.rboxlength = float(f.read(16))
        rewind += 16
        self.zboxlength = float(f.read(16))
        rewind += 16
        self.R0 = float(f.read(16))
        rewind += 16
        self.rboxleft = float(f.read(16))
        rewind += 16
        self.zbox_mid = float(f.read(16))
        f.readline()
      except:
        f.seek(-rewind, 1)
        line = f.readline().split()
        self.rboxlength = float(line[0])
        self.zboxlength = float(line[1])
        self.R0 = float(line[2])
        self.rboxleft = float(line[3])
        self.zbox_mid = float(line[4])

      self.R = ( self.rboxleft
                      + np.linspace(0, self.rboxlength, self.nrgr) )
      self.Z = ( self.zbox_mid - 0.5*self.zboxlength
                      + np.linspace(0, self.zboxlength, self.nzgr) )

      # Third line
      rewind = 0
      try:
        # Read magnetic axis
        rewind += 16
        self.Rpsi0 = float(f.read(16))  # assumed to be the same as RMAXIS
        rewind += 16
        self.Zpsi0 = float(f.read(16))  # assumed to be the same as ZMAXIS
        # Read psiaxis and psiedge
        rewind += 16
        self.PsiaxisVs = float(f.read(16))
        rewind += 16
        self.PsiedgeVs = float(f.read(16))
        # Read B0 and current
        rewind += 16
        self.Btor_at_R0 = float(f.read(16))
        f.readline()
      except:
        f.seek(-rewind, 1)
        line = f.readline().split()
        self.Rpsi0 = float(line[0])
        self.Zpsi0 = float(line[1])
        self.PsiaxisVs = float(line[2])
        self.PsiedgeVs = float(line[3])
        self.Btor_at_R0 = float(line[4])

      # Fourth line
      self.Ip = float(f.read(16))
      f.readline()

      # Fifth line
      f.readline()

      # From sixth line:

      # Read T(psi) - f (size nrgr)
      self.psieqprof = [1.0/(self.nrgr-1)*k for k in range(self.nrgr)]
      self.fprof = np.empty(self.nrgr)
      for k in range(self.nrgr):
        self.fprof[k] = readblock(f)


      # Read p(psi) (size nrgr)
      self.ptotprof = np.empty(self.nrgr)
      for k in range(self.nrgr):
        self.ptotprof[k] = readblock(f)


      # Read TT'(psi) - ff' (size nrgr)
      self.fdfdpsiprof = np.empty(self.nrgr)
      for k in range(self.nrgr):
        self.fdfdpsiprof[k] = readblock(f)


      # Read p'  (size nrgr)
      self.dpressdpsiprof = np.empty(self.nrgr)
      for k in range(self.nrgr):
        self.dpressdpsiprof[k] = readblock(f)


      # Read psi (size nrgr x nzgr in file?, nzgr x nrgr in 'output')
      self.PsiVs = np.empty(self.nrgr*self.nzgr)
      for k in range(self.nrgr*self.nzgr):
          self.PsiVs[k] = readblock(f)

      self.psinorm = [(PsiV-self.PsiaxisVs+self.PsiedgeVs)/(-self.PsiaxisVs+self.PsiedgeVs) for PsiV in self.PsiVs]
      self.PsiVs = self.PsiVs.reshape(self.nzgr, self.nrgr)

      # Read Q profile
      self.qprof = np.empty(self.nrgr)
      # self.rho_poloidal = np.empty(self.nrgr)
      for k in range(self.nrgr):
        self.qprof[k] = readblock(f)
        # self.rho_poloidal[k] = float(k)/float(self.nrgr-1)


      line = f.readline().split()
      if len(line) < 1:
        line = f.readline().split()  # CHEASE case
      if len(line) < 1:
        line = f.readline().split()  # PROCESS case
      self.npbound = int(line[0])
      self.nplimiter = int(line[1])

      self.Lcfs = np.empty(2*self.npbound)
      for k in range(2*self.npbound):
        self.Lcfs[k] = readblock(f)

      self.Lcfs = self.Lcfs.reshape(self.npbound, 2)
      self.Lcfs = self.Lcfs.transpose()

      self.Limiter = np.empty(2*self.nplimiter)
      for k in range(2*self.nplimiter):
        self.Limiter[k] = readblock(f)

      self.Limiter = self.Limiter.reshape(self.nplimiter, 2)
      self.Limiter = self.Limiter.transpose()

      try:
        # Read number of coils
        line = f.readline().split()
        if len(line) < 1:
          line = f.readline().split()  # PROCESS case
        self.ncoils = int(line[0])
      except ValueError:
        return  # No coil data available
      except IndexError:
        return  # No coil data available

      # Read PFcoilData
      self.PFcoilData = np.empty(5*self.ncoils)
      for k in range(5*self.ncoils):
        self.PFcoilData[k] = readblock(f)

      self.PFcoilData = self.PFcoilData.reshape(self.ncoils, 5)

      self.Rcond = self.PFcoilData[:,0]
      self.Zcond = self.PFcoilData[:,1]
      self.DRcond = self.PFcoilData[:,2]
      self.DZcond = self.PFcoilData[:,3]
      self.Icond = self.PFcoilData[:,4]

      try:
        # Read Brn (optional)
        self.Brn = np.empty(self.nrgr*self.nzgr)
        for k in range(self.nrgr*self.nzgr):
          self.Brn[k] = readblock(f)
        self.Brn = self.Brn.reshape(self.nzgr, self.nrgr)

        # Read Bzn (optional)
        self.Bzn = np.empty(self.nrgr*self.nzgr)
        for k in range(self.nrgr*self.nzgr):
          self.Bzn[k] = readblock(f)
        self.Bzn = self.Bzn.reshape(self.nzgr, self.nrgr)
      except (ValueError, IndexError):
        pass  # Brn/Bzn data not available


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
    from math import cos, pi, sin

    self.coilgroups = []

    self.number_toroidal_points = 128

    for k in range(len(indices_coilgroupstarts)):
      i = indices_coilgroupstarts[k]
      imax = len(self.Rcond)
      if k < len(indices_coilgroupstarts)-1:
        imax = indices_coilgroupstarts[k+1]

      for j in range(i, imax):
        c = Coilgroup(self.Icond[j], coilgroups[k], coiltags[k])

        if c.Tag in {'CS','PF'}:
          for l in range(self.number_toroidal_points):
            alpha = float(l) * 2.0*pi/float(self.number_toroidal_points)
            p = Point(self.Rcond[j]*cos(alpha),
                      self.Rcond[j]*sin(alpha),
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


  def add_tf_coilgroup(self, filename: str, current: float, coilgroup: int):
    """

    input:
    ------
    filename: string, name and path of file where to find the data.
    current: float, the current through the coil.
    coilgroup: integer, the id to use for the coil.

    output:
    -------
    None

    sideeffects:
    ------------
    Adds coilgroup to this object.
    """

    with open(filename) as f:
      l = f.readlines()

    l = [ls.strip() for ls in l]

    c = Coilgroup(current, coilgroup, "TF")
    for ls in l:
      [x, z] = ls.split()
      p = Point(float(x),
                0.0,
                float(z))
      c.coordinates.append(p)

    self.coilgroups.append(c)


  def add_rmp_coilgroup(self, filename: str, current: float, coilgroup: int):
    """

    input:
    ------
    filename: string, name and path of file where to find the data.
    current: float, the current through the coil.
    coilgroup: integer, the id to use for the coil.

    output:
    -------
    None

    sideeffects:
    ------------
    Adds coilgroup to this object.
    """

    with open(filename) as f:
      l = f.readlines()

    l = [ls.strip() for ls in l]

    c = Coilgroup(current, coilgroup, "RMP")
    for ls in l:
      [x, y, z] = ls.split()
      p = Point(float(x),
                float(y),
                float(z))
      c.coordinates.append(p)

    self.coilgroups.append(c)


  def get_used_coilgroup_ids(self):
    """
    input:
    ------
    None

    output:
    -------
    Set with all used coilgroup ids.

    sideeffects:
    ------------
    None
    """
    s = set()
    for c in self.coilgroups:
      s.add(c.coilgroup)

    return s


  def get_free_coilgroup_id(self):
    """
    input:
    ------
    None

    output:
    -------
    Free coilgroup id. Currently determined as maximum used id plus one.

    sideeffects:
    ------------
    None
    """
    s = self.get_used_coilgroup_ids()

    return (max(s)+1)


  def increase_periodicity_in_coilgroups(self, periods:int):
    r"""
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
          pass # Nothing to do here, is assumed to be already symmetric with enough points.
        else:
          if c.Tag in {'RMP'}:
            cn = Coilgroup(c.I, self.get_free_coilgroup_id(), c.Tag)
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
          f.write("{:+13.7e} {:+13.7e} {:+13.7e} {:+13.7e}\n".format(p.x, p.y, p.z, c.I))

        p = c.coordinates[0]
        f.write("{:+13.7e} {:+13.7e} {:+13.7e} {:+13.7e} {} {}\n".format(p.x, p.y, p.z, 0.0, c.coilgroup, c.Tag))

      f.write("end\n")


class Coilgroup:
  """
  Collect data for a coil(group):
  - current
  - id
  - Tag
  - coordinates
  """
  def __init__(self, I: float, coilgroup: int, Tag: str):
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


def readblock(file):
  """
  Reads a block of 16 characters from a file. Corrections in case of
  newlines and double spaces are applied.
  """
  dat = file.read(16)
  if(b'\n\n' in dat):
    file.seek(-16, 1)
    dat = file.read(18).replace(b'\n', b'')
  if(b'\n' in dat):
    file.seek(-16, 1)
    dat = file.read(17).replace(b'\n', b'')
  if(b'  ' in dat):
    file.seek(-16, 1)
    dat = file.read(17).replace(b'  ', b' ')

  return dat
