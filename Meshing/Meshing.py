class Mesh_1D:
  """
  Creates a one dimensional mesh object.
  """
  def __init__(self, nz: int, L: float):
    self.nz = nz
    self.L = L
    self.dz = L/nz
