from Meshing.Meshing import *
from Fields.Fields import *
import matplotlib.pyplot as plt

"""
Kernels:
NOTES:
(1)    most kernels indexing corresponds to "cids" on the mesh. e.g. index 0 will be cell id 0
       Otherwise we would need a map that maps indexes in matrices to cell ids which we just arent
       going to deal with in this simple 1D solver.
(2)    velocity indexing should be stored on faces of the mesh - e.g. faces between cells.
"""

### KERNEL BASE DEFINITION ###
class Kernel:
  def __init__(self, field: Field, mesh: Mesh_1D):
    self.field = field
    self.mesh = mesh
    self.b = None
    self.aC = None # matrix aF
    self.aF = None # matrix aC

    # initialize zeros --- really just initializing size of matrices
    self.b = np.zeros(self.mesh.nz)

    for idx in range(self.mesh.nz):
      try:
        self.aC = np.vstack([self.aC, np.zeros(self.mesh.nz)])
        self.aF = np.vstack([self.aF, np.zeros(self.mesh.nz)])
      except:
        self.aC = np.zeros(self.mesh.nz)
        self.aF = np.zeros(self.mesh.nz)

  def update_coeffs(self):
    self.get_aC()
    self.get_aF()
    self.get_b()
  def get_aC(self):
    pass
  def get_aF(self):
    pass
  def get_b(self):
    pass
  def plot_matrix(self):
    self.update_coeffs()
    A = abs(self.aC + self.aF)
    mask = A < 1e-10
    plt.imshow(mask, cmap='gray', vmin=0, vmax=1)
    # plt.axis('off')  # Optional: hides axis ticks
    plt.show()

  def plot_b(self):
    self.update_coeffs()
    plt.plot(self.b, 'ks-', markerfacecolor='w')



### DIFFERENT KERNELS TO USE ###
class MomentumEquation(Kernel):
  """
    A advection equation for momentum - equation 15.80 in my textbook.
  """
  def __init__(self, p: Field, mesh: Mesh_1D, w: Field, scheme: str, rho: float, mu: float):
    self.w = w
    self.p = p
    self.mesh = mesh
    self.rho = rho

    # setup for advection scheme
    self.scheme = scheme
    valid_schemes = {'upwind'}
    if self.scheme not in valid_schemes:
      raise ValueError("Scheme for advection not in allowed list of schemes")

  def get_aC(self):
    # First interpolate velocity field to faces:
    face_vals = self.w.interpolate_to_faces() # index 0 of vals corresponds to lower face of cell index 0

    # Time term Eq. 15.73
    flux_CC = 0.0 #

    # Convection + diffusion terms Eq. 15.72
    flux_Cf =

    #
  def get_aF(self):
    pass
  def get_b(self):
    pass
  def RhieChowInterpolation():
    """
      Uses rhie chow interpolation to get face values...
    """




class AdvectionKernel(Kernel):
  """
    Variable velocity advection Kernel - some scalar variable being advected by a velocity FaceField - w.
  """
  def __init__(self, field, mesh, w: FaceField, scheme: str, rho: float):
    """
    field: scalar field that is being advected.
    mesh: mesh that we are working on.
    w: z velocity we are transporting.
    scheme: which advection scheme to use. (upwind, )
    rho: fluid density (normally just 1.0 for basic advection)

    Notes:
      (1) need to make rho a cell-centered field
      (2) The velocity is stored on the faces.
    """
    # super
    super().__init__(field, mesh)

    # setup variables and values
    self.w = w
    self.rho = rho

    # setup for advection scheme
    self.scheme = scheme
    valid_schemes = {'upwind', 'quick'}
    if self.scheme not in valid_schemes:
      raise ValueError("Scheme for advection not in allowed list of schemes")


  def get_aC(self):
    # Reset aC
    self.aC *= 0.0

    ####### UPWIND SCHEME #######
    # aC = -(aE + aW) = ||-m_e, 0|| + ||-m_w, 0|| =
    if self.scheme == 'upwind':
      for cid in self.mesh.cidList:
        if self.mesh.cells[cid].upperType == 'f': # if not boundary
          aU = self.get_aUpperLowerUpwind(cid=cid, upper_or_lower='upper')
          self.aC[cid,cid] -= aU
        if self.mesh.cells[cid].lowerType == 'f':
          aL = self.get_aUpperLowerUpwind(cid=cid, upper_or_lower='lower')
          self.aC[cid,cid] -= aL
    #############################

    ######## QUICK SCHEME #######
    if self.scheme == 'quick':
      for cid in self.mesh.cidList:
        mw = self.get_m_values(cid=cid, upper_or_lower='lower')
        me = self.get_m_values(cid=cid, upper_or_lower='upper')
        # Get coeffs for this CID
        aE = -3/4*max(-me, 0) + 3/8 * max(me, 0) - 1/8*max(mw,0)
        aW = -3/4*max(-mw, 0) + 3/8 * max(mw, 0) - 1/8*max(me,0)
        aEE = 1/8*max(-me, 0)
        aWW = 1/8*max(-mw, 0)
        if self.mesh.cells[cid].upperType == 'f':
          # subtract upper cell coeff.
          self.aC[cid,cid] -= aE
          upperCid = self.mesh.cells[cid].upperNeighborCid
          if self.mesh.cells[upperCid].upperType == 'f':
            self.aC[cid,cid] -= aEE

        if self.mesh.cells[cid].lowerType == 'f':
          # subtract lower cell coeff.
          self.aC[cid,cid] -= aW
          lowerCid = self.mesh.cells[cid].lowerNeighborCid
          if self.mesh.cells[lowerCid].lowerType == 'f':
            self.aC[cid,cid] -= aWW
    #############################


  def get_aF(self):
    self.aF *= 0.0
    ####### UPWIND SCHEME #######
    if self.scheme == 'upwind':
      for cid in self.mesh.cidList:
        if self.mesh.cells[cid].upperType == 'f': # if not boundary
          aU = self.get_aUpperLowerUpwind(cid=cid, upper_or_lower='upper')
          self.aF[cid, cid+1] = aU
        if self.mesh.cells[cid].lowerType == 'f':
          aL = self.get_aUpperLowerUpwind(cid=cid, upper_or_lower='lower')
          self.aF[cid, cid-1] = aL
    #############################

    ######## QUICK SCHEME #######
    if self.scheme == 'quick':
      for cid in self.mesh.cidList:
        mw = self.get_m_values(cid=cid, upper_or_lower='lower')
        me = self.get_m_values(cid=cid, upper_or_lower='upper')
        # Get coeffs for this CID
        aE = -3/4*max(-me, 0) + 3/8 * max(me, 0) - 1/8*max(mw,0)
        aW = -3/4*max(-mw, 0) + 3/8 * max(mw, 0) - 1/8*max(me,0)
        aEE = 1/8*max(-me, 0)
        aWW = 1/8*max(-mw, 0)
        if self.mesh.cells[cid].upperType == 'f':
          # subtract upper cell coeff.
          self.aC[cid,cid+1] = aE
          upperCid = self.mesh.cells[cid].upperNeighborCid
          if self.mesh.cells[upperCid].upperType == 'f':
            self.aC[cid,cid+2] = aEE

        if self.mesh.cells[cid].lowerType == 'f':
          # subtract lower cell coeff.
          self.aC[cid,cid-1] = aW
          lowerCid = self.mesh.cells[cid].lowerNeighborCid
          if self.mesh.cells[lowerCid].lowerType == 'f':
            self.aC[cid,cid-2] = aWW
    #############################


  def get_b(self):
    pass

  def get_m_values(self, cid: int, upper_or_lower: str):
    if upper_or_lower == 'upper':
      return self.mesh.cells[cid].upperArea * self.w.get_upper(cid=cid) * self.rho # S in positive direction
    elif upper_or_lower == 'lower':
      return -1.0 * self.mesh.cells[cid].lowerArea * self.w.get_lower(cid=cid) * self.rho # S in negative direction
    else:
      raise Exception("upper or lower must be either upper or lower strings and not anything else!")

  def get_aUpperLowerUpwind(self, cid: int, upper_or_lower: str):
    """
    Returns aUpper or aLower for upwind scheme
    """
    if upper_or_lower == 'upper':
      M_upper = self.get_m_values(cid=cid, upper_or_lower='upper')
      a = -1.0 * max(-1.0*M_upper, 0.0)
    elif upper_or_lower == 'lower':
      M_lower = self.get_m_values(cid=cid, upper_or_lower='lower')
      a = -1.0 * max(-1.0*M_lower, 0.0)
    else:
      raise Exception("Must be upper or lower!")
    return a


class DiffusionKernel(Kernel):
  def __init__(self, field, mesh, Gamma: float):
    super().__init__(field, mesh)
    self.Gamma = Gamma # diffusion coefficient

  def get_aC(self):
    # assigns diagonal coefficients
    # reset aC
    self.aC *= 0.0

    # iterate for every cell now
    for cid in self.mesh.cidList:
      if self.mesh.cells[cid].upperType == 'f': # if not boundary
        self.aC[cid, cid] -= -1.0*(self.Gamma * self.mesh.cells[cid].upperNeighborGdiff)

      if self.mesh.cells[cid].lowerType == 'f': # if not boundary
        self.aC[cid, cid] -= -1.0*(self.Gamma * self.mesh.cells[cid].lowerNeighborGdiff)

  def get_aF(self):
    # assigns offdiagonal coefficients
    for cid in self.mesh.cidList:

      # assign upper
      if self.mesh.cells[cid].upperType == 'f':
        self.aF[cid, self.mesh.cells[cid].upperNeighborCid] = -self.Gamma * self.mesh.cells[cid].upperNeighborGdiff

      # assign lower
      if self.mesh.cells[cid].lowerType == 'f':
        self.aF[cid, self.mesh.cells[cid].lowerNeighborCid] = -self.Gamma * self.mesh.cells[cid].lowerNeighborGdiff

  def get_b(self):
    pass ## do nothing for this.

class ExplicitSource(Kernel):
  """
  Explicit source. Can input as a constant float, a numpy array, or as a field. If adding as a field,
  the value of T from the field will be used in an explicit fashion
  """
  def __init__(self, field: Field, mesh: Mesh_1D, Q: (np.ndarray, float, Field) ):
    super().__init__(field, mesh)

    self.Q = np.ones(len(self.b)) # sets Q to all ones

    # if it is just a float --
    if isinstance(Q, float):
      self.Q *= Q
    # if it is a numpy array, just set equal to the reference of that array
    elif isinstance(Q, np.ndarray):
      self.Q = Q

    # if it is a field, now set self.Q equal to the field so when we reference it later, self.Q will be Field
    elif isinstance(Q, Field):
      self.Q = Q
    else:
      raise ValueError("Q type must be np array, float, or a Field for explicit source kernel")

  def get_b(self):
    self.b *= 0.0
    if isinstance(self.Q, Field):
      self.b += self.Q.T # gets value from the field
    else:
      self.b += self.Q # else it just adds it like a np array or float array.



### BOUNDARY CONDITIONS ###
class BoundaryCondition(Kernel):
  def __init__(self, field, mesh, boundary: str):
    super().__init__(field, mesh)

    # basic checks
    if (boundary is not 'upper') & (boundary is not 'lower'):
      raise Exception("Boundary is not upper or lower")

    # assign upper or lower
    self.boundary = boundary

class DirchletBC(BoundaryCondition):
  def __init__(self, field, mesh, boundary, phi: float, Gamma: float):
    super().__init__(field, mesh, boundary)

    # assignments
    self.phi = phi
    self.Gamma = Gamma

  def get_face_value(self):
    """
      gets value at the associated face to use when computing phi_f
      used for interpolation - especially with gradients and constructing interpolated face values
    """
    return self.phi

  def get_face_gradient(self):
    """
      gets gradient value at the face for this boundary - useful for gradient stuffs.
    """
    return 0.0

  def get_b(self):
    if self.boundary == 'upper': # upper boundary condition
      cid = self.mesh.cidList[-1] # last
      Sb = self.mesh.cells[cid].upperArea

    if self.boundary == 'lower': # first cell
      cid = self.mesh.cidList[0]
      Sb = self.mesh.cells[cid].lowerArea

    dCb = self.mesh.cells[cid].dz / 2.0
    self.b[cid] = 1.0 * self.Gamma * Sb / dCb * self.phi # b = -FluxVb where fluxVb = -Gamma_b * gDiff_b = -Gamma_b * S_b / dCb

  def get_aC(self):
    if self.boundary == 'upper': # upper boundary condition
      cid = self.mesh.cidList[-1] # last
      Sb = self.mesh.cells[cid].upperArea
    if self.boundary == 'lower': # first cell
      cid = self.mesh.cidList[0]
      Sb = self.mesh.cells[cid].lowerArea

    dCb = self.mesh.cells[cid].dz / 2.0
    self.aC[cid, cid] = 1.0 * self.Gamma * Sb / dCb # b = -FluxVb where fluxVb = -Gamma_b * gDiff_b = -Gamma_b * S_b / dCb

class InletBC(BoundaryCondition):
  """
  Boundary condition for the inlet of the problem.
  """
  def __init__(self, field, mesh, boundary: str):
    super().__init__(field, mesh, boundary)
  def get_face_gradient(self):
    pass
  def get_face_value(self):
    pass
  def get_aC(self):
    pass
  def get_aF(self):
    pass
  def get_b(self):
    pass



